import pandas as pd


class ResultsQueries(object):
	def __init__(self,results):
		self._results = results
		duration = self._results.runDetails.end - self._results.runDetails.startRecording
		self._years = duration.days/365.25
		self._convert = {}
		self._convert['kg'] = 1
		self._convert['t'] = 0.001
		self._convert['y'] = 1/self._years
		self._convert['t/y'] = self._convert['t'] * self._convert['y']
		self._convert['kg/y'] = self._convert['kg'] * self._convert['y']
		self._convert['standard'] = {
			'Sediment - Fine':1e-6,'Flow':1e-6,'N_DIN':1e-3,'N_DON':1e-3,
			'N_Particulate':1e-3,'P_DOP':1e-3,'P_FRP':1e-3,'P_Particulate':1e-3
		}

	def convert(self,data,units):
		conversion = self._convert[units]
		if hasattr(conversion,'keys'):
			result = data.copy()
			for k,v in conversion.items():
				result[k] *= v
			return result
		else:
			return data * conversion

	def regional_export(self,units='t/y'):
		full = self._results.get('RegionalSummaryTable')
		exports = full[full.MassBalanceElement=='Export']
		regionalLoads = exports.reset_index().pivot('Constituent','SummaryRegion','Total_Load_in_Kg')
		return self.convert(regionalLoads,units)

	def export_per_landuse_per_ha(self,load_units='t/y'):
		fuData_kg = self._results.get('FUSummaryTable')
		SupplyByFU_t = self.convert(fuData_kg.reset_index().pivot('Constituent','FU','Total_Load_in_Kg'),load_units)

		FUAreas_ha = self._results.get('fuAreasTable').groupby('FU').sum() / 10000.0

		return SupplyByFU_t / FUAreas_ha.Area

	def annual_mass_balance(self,units='standard'):
		raw = self._results.get('RawResults')
		grouped_rawResults = raw.reset_index().groupby(['Process', 'BudgetElement', 'Constituent']).sum()
		mass_balance_kg = pd.pivot_table(grouped_rawResults.reset_index(),index = ['Process','BudgetElement'], columns = 'Constituent', values = 'Total_Load_in_Kg')

		mass_balance_annual = self.convert(mass_balance_kg,units) / self._years

		col_order = ['Sediment - Fine','Sediment - Coarse','P_Particulate','P_FRP','P_DOP','N_Particulate','N_DIN','N_DON']
		other_cols = list(set(mass_balance_annual.columns)-set(col_order))
		mass_balance_annual = mass_balance_annual[col_order + other_cols]

		transposed=mass_balance_annual.transpose()
		elements = ['Supply','Loss','Residual','Other']
		for element in elements:
			transposed['Total',element] = mass_balance_annual.ix[element].fillna(0).sum()

		transposed['Total','Export'] = 2*transposed['Total','Supply']-transposed['Total'].transpose().sum()

		# Then transpose back...
		transposed = transposed.transpose()
		return transposed.ix[['Supply','Loss','Other','Residual','Total']]

	def compare_mean_annual_loads(self,observed_fn,monitoring_sites):
		import pandas as pd
		raw = self._results.get('RawResults')
		linkYields = raw[raw.BudgetElement=='Link Yield']
		relevantLinkYields = linkYields[linkYields.ModelElement.isin(monitoring_sites)]
		loadsForComparison = relevantLinkYields.reset_index().pivot('ModelElement','Constituent','Total_Load_in_Kg')
		conversions = {'Sediment - Fine':1e-6,'Flow':1e-6,'N_DIN':1e-3,'N_DON':1e-3,
			'N_Particulate':1e-3,'P_DOP':1e-3,'P_FRP':1e-3,'P_Particulate':1e-3}
		for col in loadsForComparison.columns:
			if col in conversions:
				loadsForComparison[col] *= conversions[col]
		
		loadsForComparison /= self._years

		observedLoads = pd.DataFrame.from_csv(observed_fn)	
		loadsForComparison['ValueType']='Predicted'
		observedLoads['ValueType']='Observed'
		return pd.concat([loadsForComparison,observedLoads]).reset_index().set_index(['ModelElement','ValueType'])

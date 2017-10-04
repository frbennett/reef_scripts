
def is_event(timeclass):
    """
    It's an event if it starts with a number...
    """
    if not isinstance(timeclass,str):
        return False
    return timeclass[0].isdigit()

def is_baseflow(timeclass):
    if not isinstance(timeclass,str):
        return False
    return timeclass.startswith('b')

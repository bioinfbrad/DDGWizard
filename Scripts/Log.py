from Scripts.Global_Value import Log_Path
def Init_Log():
    with open(Log_Path, 'w') as log:
        pass


def Log(text):
    if not isinstance(text,str):
        return
    with open(Log_Path,'a') as log:
        log.write(text+'\n')

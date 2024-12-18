from Scripts.Log import Log
class Error_Class:
    # def Chain_Num_Not_Equal_1(self, function_name,pdb):
    #     print(function_name+' - '+pdb+' - '+'Chain Number is not equal with 1')

    def Is_Not_Existed(self,function_name, what):
        Log(function_name + ' - ' + what + ' - ' + 'Which is not existed')
        print(function_name + ' - ' + what + ' - ' + 'Which is not existed')

    def Request_Error(self, function_name, what):
        Log(function_name+' - '+what+' - '+'Request of that returns error')
        print(function_name+' - '+what+' - '+'Request of that returns error')

    def Modelling_Fail(self, function_name, what):
        Log(function_name + ' - ' + what + ' - ' + 'Modelling in this function is failed')
        print(function_name + ' - ' + what + ' - ' + 'Modelling in this function is failed')

    def Something_Wrong(self,function_name, what='foo'):
        Log(function_name + ' - ' + what + ' - ' + 'In this function, something is wrong')
        print(function_name + ' - ' + what + ' - ' + 'In this function, something is wrong')


error_obj=Error_Class()
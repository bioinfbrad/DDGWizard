import dill
import pandas as pd
import joblib
import xlwt
import Scripts.Global_Value
from Scripts.Log import Log


def load_pkl(filepath):  # load model of pkl format from filepath, and return data (model)
    with open(filepath, "rb") as fr:
        data = dill.load(fr, encoding="utf-8")
    Log(f"[{filepath}] data loading...")
    print(f"[{filepath}] data loading...")
    return data


def XGBoostRegression_Predict(csv_path,model_path,output_path):
    data = pd.read_csv(csv_path)  # 读取数据
    ids=list(data['ID'])
    data = data.drop("ID", axis=1)  # 删除ID列

    rfe_infos = pd.read_excel(model_path+'rfe_infos.xlsx')
    X_cols = rfe_infos[rfe_infos["ranking"] == 1]["feature_names"].tolist()
    fea=data[X_cols].values
    sc = load_pkl(model_path+'sc.pkl')
    fea = sc.transform(fea)
    model = load_pkl(model_path+'predictor.pkl')
    y_pred = model.predict(fea)
    res_dict={}
    if not Scripts.Global_Value.Is_Use_Reverse_Data:
        for i in range(len(y_pred)):
            res_dict[ids[i]]=['forward',y_pred[i]]

        Log(res_dict)
        print(res_dict)
        Log(f'result saving in {output_path}')
        print(f'result saving in {output_path}')
        wb=xlwt.Workbook()
        ws=wb.add_sheet('sheet1')
        header=['ID','Forward/Reverse','ddG']
        for i in range(len(header)):
            ws.write(0,i,header[i])
        count=1
        for name in ids:
            tag=res_dict[name][0]
            pred_ddG=float(res_dict[name][1])
            ws.write(count,0,name)
            ws.write(count,1,tag)
            ws.write(count,2,pred_ddG)
            count+=1
        wb.save(f'{output_path}Pred_ddG.xls')
    else:
        assert len(y_pred) == len(ids)
        is_pass_rev=False
        for i in range(len(y_pred)):
            if i%2==0:
                if i+1 < len(y_pred):
                    for_id=str(ids[i])
                    for_mut=for_id.split('_')[2]
                    rev_id=str(ids[i+1])
                    rev_mut=rev_id.split('_')[2]
                    for_mut_changed=list(for_mut)
                    prefix=for_mut_changed.pop(0)
                    postfix=for_mut_changed.pop(-1)
                    mid=''.join(for_mut_changed)
                    for_mut_changed_str=postfix+mid+prefix
                    if for_id.split('_')[0]==rev_id.split('_')[0] and for_id.split('_')[1]==rev_id.split('_')[1] and for_id.split('_')[3]==rev_id.split('_')[3] and rev_mut==for_mut_changed_str:
                        res_dict[ids[i]] = ['forward',y_pred[i]]
                    else:
                        res_dict[ids[i]] = ['unknown', y_pred[i]]
                        is_pass_rev=True
                else:
                    res_dict[ids[i]] = ['unknown', y_pred[i]]
                    is_pass_rev=True
            else:
                if is_pass_rev:
                    is_pass_rev=False
                else:
                    res_dict[ids[i]] = ['reverse', y_pred[i]]

        Log(res_dict)
        print(res_dict)
        Log(f'result saving in {output_path}')
        print(f'result saving in {output_path}')
        wb=xlwt.Workbook()
        ws=wb.add_sheet('sheet1')
        header=['ID','Forward/Reverse','ddG']
        for i in range(len(header)):
            ws.write(0,i,header[i])
        count=1
        for name in ids:
            tag=res_dict[name][0]
            pred_ddG=float(res_dict[name][1])
            ws.write(count,0,name)
            ws.write(count,1,tag)
            ws.write(count,2,pred_ddG)
            count+=1
        wb.save(f'{output_path}Pred_ddG.xls')


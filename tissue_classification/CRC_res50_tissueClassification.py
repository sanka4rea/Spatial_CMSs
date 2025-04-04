import pandas as pd
from torch.utils.data import DataLoader
import torch
import torch.nn as nn
from torchvision import datasets, models, transforms
import os
# from numba import jit
import csv

BATCH_SIZE = 800
NUM_classes = 9
criterion = nn.CrossEntropyLoss()
Result_PATH = '/data/classification/'
os.environ['CUDA_VISIBLE_DEVICES'] = '4'
# torch.cuda.set_device(0, 1, 2, 3)

#数据集函数
data_transforms = {
    'val': transforms.Compose([
        transforms.Resize((224, 224)),
        transforms.ToTensor()
        #transforms.Normalize([0.485, 0.456, 0.406], [0.229, 0.224, 0.225])
    ]),
}

val_dataset = datasets.ImageFolder(root='/data/ARGO/div_220/',
                transform=data_transforms['val'])
val_D = DataLoader(val_dataset, batch_size=BATCH_SIZE, shuffle=False, num_workers=10)

device = torch.device("cuda")
model_ft = models.resnet50(pretrained=True) # , progress=True
num_ftrs = model_ft.fc.in_features
model_ft.fc = nn.Linear(num_ftrs, NUM_classes)
model_ft.load_state_dict(torch.load('./CRC_res50_tissueClass_weights_24.pt'))
model_ft = nn.DataParallel(model_ft)
net = model_ft.to(device)

net.eval()
correct = 0
total = 0
val_loss = 0.0

# 
with torch.no_grad():
    for i, data in enumerate(val_D, 0):
        images, labels_val = data
        images, labels_val = images.to(device), labels_val.to(device)
        outputs = net(images)
        # loss = criterion(outputs, labels_val)
        # val_loss += loss.item()

        _, predicted = torch.max(outputs.data, 1)
        # total += labels_val.size(0)
        # correct += (predicted == labels_val).sum().item()

        prob = outputs.data.cpu()
        pre_label = predicted.cpu()
        tr_label = labels_val.cpu()
        if i == 0:
            df = pd.DataFrame({'Val predict class': list(pre_label.numpy()), 'True label': list(tr_label.numpy())}) # , index=[0]
            # df_prob = pd.DataFrame(data=prob.numpy(), index=range(len(prob)))
        else:
            df_tmp = pd.DataFrame({'Val predict class': list(pre_label.numpy()), 'True label': list(tr_label.numpy())})
            # df_prob_tmp = pd.DataFrame(data=prob.numpy(), index=range(len(prob)))
            df = pd.concat([df, df_tmp])
            # df_prob = pd.concat([df_prob, df_prob_tmp], ignore_index=True)

        # print("Finish ", i, "Batches")

df.to_csv(os.path.join(Result_PATH, "predict.csv"), index=True, sep=',')
test=pd.DataFrame(data=val_dataset.imgs)
test.to_csv(os.path.join(Result_PATH, "True.csv"), index=True, sep=',')
import pandas as pd

# 读取CSV文件，不使用表头
df = pd.read_csv('spikeCount.csv', header=None)
print(df.shape)

# 统计每1000行中大于15的数据个数
count = 0
for i in range(0, len(df), 1000):
    subset = df[i:i+1000]  # 获取每1000行的子集
    count += (subset > 10).sum().sum()  # 统计大于15的数据个数

print("大于15的数据个数：", count)

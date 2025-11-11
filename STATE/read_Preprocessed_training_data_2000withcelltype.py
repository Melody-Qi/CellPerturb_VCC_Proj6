import h5py

f = h5py.File("../vcc_data/preprocessed_training_data_2000withcelltype.h5ad", "r")
print(list(f["obs"].keys()))  # 看看 obs 里有哪些键

# 打印 cell_type 的内容
print(f["obs/cell_type/categories"][:])  # 类别名
print(f["obs/cell_type/codes"][:10])     # 前10个 cell 的类别编号

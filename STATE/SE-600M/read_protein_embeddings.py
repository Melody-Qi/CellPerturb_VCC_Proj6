import torch

path = "protein_embeddings.pt"   # 修改成你的文件路径

print(f"Loading: {path}")
obj = torch.load(path, map_location="cpu")

print("\n=== Object Type ===")
print(type(obj))

if isinstance(obj, dict):
    print("\n=== Keys in dict ===")
    print(list(obj.keys())[:20])  # 只打印前20个键

    # 尝试打印每个 key 对应对象的尺寸
    for k, v in list(obj.items())[:5]:
        print(f"\nKey: {k}")
        if torch.is_tensor(v):
            print(f" Tensor shape: {v.shape}")
        else:
            print(f" Value type: {type(v)}")

elif torch.is_tensor(obj):
    print("\n=== Tensor details ===")
    print("shape:", obj.shape)
    print("dtype:", obj.dtype)

elif isinstance(obj, list):
    print("\n=== List details ===")
    print("length:", len(obj))
    if len(obj) > 0 and torch.is_tensor(obj[0]):
        print("First element shape:", obj[0].shape)

else:
    print("\nUnknown Type Content:")
    print(obj)

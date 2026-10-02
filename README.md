# task26100305

含重复值的 Mann-Whitney 精确检验与位置差（A − B）置信集合反演库。
输入为非空、已清洗的有限十进制字符串序列，一律按精确有理数解释，
全部计算基于 `fractions.Fraction`，端点与 p 值完全精确。

## 调用

```python
from mwexact import exact_test, invert_shift

# 精确检验：返回观测统计量 U 与精确双侧 p 值（均为 Fraction）
result = exact_test(["1", "2", "4", "5"], ["0", "1", "3"])
print(result.u, result.p)          # U 与精确 p 值

# 置信集合反演：置信度支持 90% / 95% / 99%
cs = invert_shift(["1", "2", "4", "5"], ["0", "1", "3"], 0.95)
print(cs)                          # 精确区间并集，如 (-2, 5)
cs.contains("1.5")                 # 任意十进制字符串 / int / Fraction
cs.intervals                       # Interval(low, high, low_closed, high_closed)
```

- 置信度可写成 `0.9`、`"90%"`、`"0.95"`、`Fraction(19, 20)` 或 `95`。
- 置信集合恰好包含 p 值严格大于显著性水平的 theta；p 恰好等于
  显著性水平的 theta 被排除（边界为开）。
- 区间端点为 `Fraction` 或 `None`（无界）；单点集合表示为
  `low == high` 且两侧均闭。

## 统计口径

- U：跨组对 A>B 计 1、相等计 1/2、小于计 0，求和。
- 原假设：合并样本重新分配、组大小不变，每种分配等可能；
  数值相同的观测仍按不同个体计数。
- 双侧 p：重分配后 |U − nm/2| 不小于观测距离的分配比例。
- 反演：对 A 中每个值减去 theta 后重新检验，保留 p > 1 − 置信度
  的全部 theta。

## 测试

```sh
python3 -m unittest discover -s tests -v
```

测试以独立的标签枚举参照（`tests/reference.py`）核对：
含重复值且组大小不等样本的 U 与精确 p；置信集合在跨组差值本身、
相邻差值中点及两侧外部的成员关系；以及 p 恰好等于显著性水平时
端点被排除的情形。

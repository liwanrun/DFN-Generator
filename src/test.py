import matplotlib.pyplot as plt
import numpy as np
from shapely import LineString
from shapely import Polygon
from fracture import RoughFracture

# 创建一个多段线对象
frac = RoughFracture()
dn, phi = frac._calc_harmonic()
frac._transform_IFFT(dn, phi)
line = LineString(frac.Points)

# 计算多段线的envelope
buffer = line.buffer(0.1)
print(type(buffer))

# 获取多段线和envelope的坐标
line_coords = np.array(line.coords)
envelope_coords = np.array(buffer.exterior.coords)

# 绘制多段线和envelope
plt.figure()
plt.plot(line_coords[:, 0], line_coords[:, 1], label='LineString')
plt.plot(envelope_coords[:, 0], envelope_coords[:, 1], label='Envelope', linestyle='--')
plt.xlabel('X')
plt.ylabel('Y')
plt.gca().set_ylim([-1.0, 1.0])
plt.legend()
plt.title('LineString and its Envelope')
plt.show()
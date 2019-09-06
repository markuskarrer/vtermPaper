#script calculating the ISA (International Standard atmosphere)

import sys
import math

g = 9.80665
R = 287.00
def cal(p0, t0, a, h0, h1):
	if a != 0:
		t1 = t0 + a * (h1 - h0)
		p1 = p0 * (t1 / t0) ** (-g / a / R)
	else:
		t1 = t0
		p1 = p0 * math.exp(-g / R / t0 * (h1 - h0))
	return t1, p1

def isa(altitude):
	a = [-0.0065, 0, 0.001, 0.0028]
	h = [11000, 20000, 32000, 47000]
	p0 = 101325
	t0 = 288.15
	prevh = 0
	if altitude < 0 or altitude > 47000:
		print("altitude must be in [0, 47000]")
		return
	for i in range(0, 4):
		if altitude <= h[i]:
			temperature, pressure = cal(p0, t0, a[i], prevh, altitude)
			break;
		else:
			# sth like dynamic programming
			t0, p0 = cal(p0, t0, a[i], prevh, h[i])
			prevh = h[i]

	density = pressure / (R * temperature)
	strformat = 'Temperature: {0:.2f} \nPressure: {1:.2f} \nDensity: {2:.4f}'
	#print(strformat.format(temperature, pressure, density))
	return(temperature, pressure, density)

if __name__ == '__main__':
	# command line option
	if len(sys.argv) != 2:
		print("Usage: ./isa.py <altitude>")
		sys.exit(-1)
	try:
		altitude = float(sys.argv[1])
	except:
		print("Wrong altitude format!")
		sys.exit(-2)
        isa(altitude)
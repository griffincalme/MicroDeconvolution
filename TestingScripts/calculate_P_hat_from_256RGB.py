import math

red = 236
green = 32
blue = 80

print('Original 8 bit channel RGB values:')
print(red)
print(green)
print(blue)
print('------------\n')

print('Decimal RGB:')
print(red/255)
print(green/255)
print(blue/255)
print('------------\n')



def OD_from_Ic(rgb_255):
    rgb_dec = rgb_255 / 255
    OD = -1 * math.log10(rgb_dec)
    return OD


p_red = OD_from_Ic(red)
p_green = OD_from_Ic(green)
p_blue = OD_from_Ic(blue)


print('Optical Density of RGBs:')
print(p_red)
print(p_green)
print(p_blue)
print('------------\n')

p_squared_sum = (p_red ** 2) + (p_green ** 2) + (p_blue ** 2)

p_hat_red = p_red / math.sqrt(p_squared_sum)
p_hat_green = p_green / math.sqrt(p_squared_sum)
p_hat_blue = p_blue / math.sqrt(p_squared_sum)


print('Normalized values:')
print(p_hat_red)
print(p_hat_green)
print(p_hat_blue)


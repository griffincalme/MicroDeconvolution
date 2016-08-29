import math
import numpy as np


#Example for Hematoxylin, Eosin, DAB
# R      G      B
# X      X      X   Hematoxylin(0)
# X      X      X   Red(1)
# X      X      X   DAB(2)


OD_matrix = np.array([[0.18, 0.20, 0.08],
                      [0.1187, 0.5432, 0.6661],
                      [0.10, 0.21, 0.29]])

normalized_OD_matrix = np.array([0,0,0])


print('original')
print(OD_matrix)
print('_________')

print('\n')
print('normalized')


def get_p_hat_row(p_red, p_green, p_blue):

    p_squared_sum = p_red ** 2 + p_green ** 2 + p_blue ** 2

    p_hat_red = p_red / math.sqrt(p_squared_sum)
    p_hat_green = p_green / math.sqrt(p_squared_sum)
    p_hat_blue = p_blue / math.sqrt(p_squared_sum)

    normalized_row = [p_hat_red, p_hat_green, p_hat_blue]
    return normalized_row


for row in OD_matrix:
    p_red = row[0]
    p_green = row[1]
    p_blue = row[2]

    normalized_row = get_p_hat_row(p_red, p_green, p_blue)
    print(normalized_row)


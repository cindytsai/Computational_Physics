import numpy as np


def func1(N):

    step = np.array([[0] * (2 * N + 1)] * (2 * N + 1))
    step[N][N] = 1

    for i in range(N):
        temp_step = np.array([[0] * (2 * N + 1)] * (2 * N + 1))
        for w in range(1, np.max(step) + 1):
            coor_list = np.argwhere(step == w)
            for j in range(len(coor_list)):
                tempX = coor_list[j][0]
                tempY = coor_list[j][1]
                temp_step[tempX + 1][tempY] = temp_step[tempX + 1][tempY] + w
                temp_step[tempX - 1][tempY] = temp_step[tempX - 1][tempY] + w
                temp_step[tempX][tempY + 1] = temp_step[tempX][tempY + 1] + w
                temp_step[tempX][tempY - 1] = temp_step[tempX][tempY - 1] + w
        step = temp_step

    prob = step
    # print("RESULT:")
    # print(prob)
    # print(np.sum(prob))

    sum_distance = 0.0
    for i in range(1, int(np.max(prob) + 1)):
        coor_list = np.argwhere(prob == i)
        for j in range(len(coor_list)):
            distance = coor_list[j] - N
            distance = np.sqrt(np.sum(distance**2))
            sum_distance = sum_distance + distance * i

    mean_distance = sum_distance / np.sum(prob)
    print (N, "   ",mean_distance)


N = 1
while (N <= 101):
    # print(N)
    func1(N)
    N = N + 1

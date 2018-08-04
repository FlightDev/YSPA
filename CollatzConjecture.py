

dict = {}

nums = []

for i in range(2, 1000):
    temp = i
    counter = 0
    while temp >= i:
        if temp % 2 == 0:
            temp /= 2
        else:
            temp = temp * 3 + 1
        counter += 1
    nums.append(counter)

print nums

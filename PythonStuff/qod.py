import matplotlib.pyplot as plt
"""
def to_decimal_from_four(s):
    sum = 0
    for i in range(len(s)):
        sum += int(s[i]) * (4 ** i)
    return sum

file = open("qod.txt", "r")
s = ""
for line in file:
    s += line
nums = s.split(" ")
for i in range(len(nums)):
    if i == len(nums) - 1:
        nums[i] = nums[i][0:len(nums)-2]
    nums[i] = '0b' + nums[i]
digs = []
for num in nums:
    digs.append(int(num, 2))
final = ""
sum = 0
print min(digs), max(digs), len(digs)
nums = []
word = ""
for i in digs:
    word += str(i-48)
    nums.append(i-48)

print nums
print word
words = word.split("3")
print words
for i in words:
    print i.split("2")

zero = []
one = []
two = []
three = []
for i in range(len(nums)):
    if nums[i] == 0:
        zero.append(i)
    elif nums[i] == 1:
        one.append(i)
    elif nums[i] == 2:
        two.append(i)
    elif nums[i] == 3:
        three.append(i)

print zero
print one
print two
print three
sum = len(zero)
for i in range(len(one)):
    sum += one[i]
for i in range(len(two)):
    sum += (two[i] ** 2)
for i in range(len(three)):
    sum += (two[i] ** 3)

print sum
print len(zero), len(one), len(two), len(three)
"""
"""
i = 1
for let in digs:
    final += str(unichr(let))
    if i % 12 == 0:
        final += "\n"
    i += 1
fours = []
temp = []
j = 1
print final
for i in range(len(final)):
    temp.append(final[i])
    if j % 4 == 0:
        fours.append(temp)
    j += 1
lines = final.split("\n")
print lines
print len(lines[0])
c0x = []
c0y = []
c1x = []
c1y = []
c2x = []
c2y = []
c3x = []
c3y = []
cx = []
cy = []
sum = 0

for line in range(len(lines) - 1):
    for c in range(len(lines[line])):
        if int(lines[line][c]) == 0:
            c0x.append(line)
            c0y.append(c)
        elif int(lines[line][c]) == 1:
            c1x.append(line)
            c1y.append(c)
        elif int(lines[line][c]) == 2:
            c2x.append(line)
            c2y.append(c)
        elif int(lines[line][c]) == 3:
            c3x.append(line)
            c3y.append(c)
        else:
            print "U DUMMY"
        sum += 28 ** int(lines[line][c])

print sum


plt.style.use('dark_background')
plt.scatter(c0x, c0y, color=(1, 0, 0, 1.0))
plt.scatter(c1x, c1y, color=(0, 1, 0, 1.0))
plt.scatter(c2x, c2y, color=(0, 0, 1, 1.0))
plt.scatter(c3x, c3y, color=(1.00, 1.00, 1.00, 1.0))
"""
def find_words(s):
    temp = []
    for i in range(26):
        temp.append((chr(i + 65) + s + chr(i + 65)))
    return temp
word = ["ROU", "ARSNI", "HRIF", "REGAN", "YPIS", "UMM", "VERD", "OTATO", "OILE", "AYA", "RUS", "ARDIA"]
w = []
for i in word:
    w.append(find_words(i))

for i in range(len(w[0])):
    s = ""
    for j in range(len(w)):
        s += w[j][i] + "\t"
    print s

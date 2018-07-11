from visual import *
from datetime import datetime

#ring
#arrow
#sphere
original_hour = vector(0, 5, 0)
original_minute = vector(0, 8, 0)
original_second = vector(0, 10, 0)
clock_face = ring(pos = (0, 0, 0), axis = (0, 0, 1), thickness = 0.5, radius = 10, color = color.white)
top_label = label(pos = (0, 10, 0), text = "12", color = color.white, height = 32, box = False)
middle_thing = sphere(pos = (0, 0, 0), radius = 1)
hour_hand = arrow(pos = (0, 0, 0), axis = original_hour, shaftwidth = 0.5, color = color.yellow)
minute_hand = arrow(pos = (0, 0, 0), axis = original_minute, shaftwidth = 0.4, color = color.green)
second_hand = arrow(pos = (0, 0, 0), axis = original_second, shaftwidth = 0.2, color = color.red)

def rotate_to_time(time):
    seconds = (time[2]) * 6.0
    minutes = (time[1]) * 6.0
    hours = (time[0]) * 30
    #print seconds, minutes, hours
    seconds = radians(seconds)
    minutes = radians(minutes)
    hours = radians(hours)
    print seconds, minutes, hours
    second_hand.axis = rotate(vector(0, 10, 0), -1 * seconds, vector(0, 0, 1))
    minute_hand.axis = rotate(vector(0, 8, 0), -1 * minutes, vector(0, 0, 1))
    hour_hand.axis = rotate(vector(0, 5, 0), -1 * hours, vector(0, 0, 1))



time = [datetime.now().hour, datetime.now().minute, int(datetime.now().second)]

while True:
    rate(1)
    time[2] += 1
    if time[2] == 60:
        time[2] = 0
        time[1] += 1
    if time[1] == 60:
        time[1] = 0
        time[0] += 1
    if time[0] == 24:
        time[0] == 0
    print time
    rotate_to_time(time)

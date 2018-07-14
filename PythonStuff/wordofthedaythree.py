"""
def fun(l):
    s = ""
    for i in range(324):
        s += l[i]
        if (i + 1) % 18 == 0:
            s += "\n"
        else:
            s += " "
    return s

print fun("DROFNARBSAEEFRUDEEBASSPGBEINECKEHMMSNTHGIELKOORBYASASLORVCIMVLLETTABSLWLOUOEHCABINETSZSOIARMBSSHLORIAIPQUNMGKBKSEROMSICAMPZECNMUFEMNOPGNYPSHVOCIHLRNAFKMPNKEIBAWLRMLAKSLQREEBALINDEOSRNROOLWAERBLNDALSTQKANVHLRNQOIGEVLEERLHSIIKDCZDMHREANROIQTNEECULYAABNNFLSNNSLTEERTSNMIPREIEELEGDIRBDOOWLOZLNYEYICFWOOLSEYTREDGYRKARHOHPLODURTUX")

solutions = ['ARTHUR', 'WATSON', 'BASS', 'BOYER', 'BECTON', 'BENJAMIN', 'FRANKLIN', 'BERKELEY', 'CHARLES', 'BINGHAM', 'BRADY', 'BRANFORD', 'BEINECKE', 'RARE', 'CONNECTICUT', 'CLASS', 'DURFEE', 'DAVIES',  'BECTON', 'DAVENPORT', 'DUNHAM', 'DOW', 'EDWIN', 'MCCLELLAN', 'EZRA', 'STILES', 'EDWARD', 'EVANS', 'GRACE', 'HOPPER', 'GREELEY', 'HOLCOMBE', 'GREEN', 'HENDRIE', 'JONATHAN', 'EDWARDS', 'KIRTLAND', 'KLINE', 'KROON', 'LAWRANCE', 'LINSLY', 'CHITTENDEN', 'LEITNER', 'ABBY', 'MITCH', 'LEIGH', 'LEET', 'OLIVER', 'JEFFREY', 'LORIA', 'HENRY', 'LUCE', 'LANMAN', 'WRIGHT', 'MORSE', 'MALONE', 'MASON', 'PAULI', 'MURRAY', 'OSBORN', 'PIERSON', 'PHELPS', 'PEABODY', 'PAYNE', WHITNEY', 'RUDOLP', 'ROSENKRANZ', 'SAGE', 'STERLING', 'SILLIMAN', 'SPRAGUE', 'SLOANE', 'SHEFFIELD', 'STRATHCONA', 'STOECKEL', 'SAYBROOK', 'ANLYAN', 'TRUMBULL', 'TIMOTHY', 'DWIGHT', 'VANDERBILT', 'WELCH', 'WRIGHT', 'WILLIAM', 'HARKNESS', 'WATSON', 'CABINET']
"""

solutions = ['ROSENKRANZ']
# Here goes your puzzle, separate each letter by a space.

puzzle ='''D R O F N A R B S A E E F R U D E E
B A S S P G B E I N E C K E H M M S
N T H G I E L K O O R B Y A S A S L
O R V C I M V L L E T T A B S L W L
O U O E H C A B I N E T S Z S O I A
R M B S S H L O R I A I P Q U N M G
K B K S E R O M S I C A M P Z E C N
M U F E M N O P G N Y P S H V O C I
H L R N A F K M P N K E I B A W L R
M L A K S L Q R E E B A L I N D E O
S R N R O O L W A E R B L N D A L S
T Q K A N V H L R N Q O I G E V L E
E R L H S I I K D C Z D M H R E A N
R O I Q T N E E C U L Y A A B N N F
L S N N S L T E E R T S N M I P R E
I E E L E G D I R B D O O W L O Z L
N Y E Y I C F W O O L S E Y T R E D
G Y R K A R H O H P L O D U R T U X
'''

# Just formats the puzzle into a more computer-readable text
wordgrid = puzzle.replace(' ','')

# Computers start counting at zero, so...
length = wordgrid.index('\n')+1


characters = [(letter, divmod(index, length))
            for  index, letter in enumerate (wordgrid)]

wordlines = {}
# These next lines just  directions so you can tell which direction the word is going
directions = {'going downwards':0, 'going downwards and left diagonally':-1, 'going downwards and right diagonally':1}

for word_direction, directions in directions.items():
    wordlines[word_direction] = []
    for x in range(length):
        for i in range(x, len(characters), length + directions):
            wordlines[word_direction].append(characters[i])
        wordlines[word_direction].append('\n')

# Nice neat way of doing reversed directions.
wordlines['going right'] = characters
wordlines['going left'] = [i for i in reversed(characters)]
wordlines['going upwards'] = [i for i in reversed(wordlines['going downwards'])]
wordlines['going upwards and left diagonally'] = [i for i in reversed(wordlines['going downwards and right diagonally'])]
wordlines['going upwards and right diagonally'] = [i for i in reversed(wordlines['going downwards and left diagonally'])]


def printitout(direction, tuple, lines):
    print "Keep in mind, rows are horizontal and columns are vertical.\n"
    for direction, tuple in lines.items():
        string = ''.join([i[0] for i in tuple])
        for word in solutions:
            if word in string:
                coordinates = tuple[string.index(word)][1]
                print word, 'is at row', coordinates[0]+1, 'and column', coordinates[1]+1, direction + "."

printitout(word_direction, tuple, wordlines)

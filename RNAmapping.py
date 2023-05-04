#!/usr/bin/env python3
def RNAmap():
	import sys
	import copy
	import random
	import os
	import shutil
	
	ex = 0
	MSL = 'default'
	MCN = 'default'
	min_loop = 'default'
	LS = 'default'
	acceptable_unmatch_num = 'default'
	
	for argv_i in range(2,7) :
	    try :
	        if sys.argv[argv_i].split(':')[0] == 'MinLength' :
	            MSL = int(sys.argv[argv_i].split(':')[1])
	        elif sys.argv[argv_i].split(':')[0] == 'MinComplement' :
	            MCN = int(sys.argv[argv_i].split(':')[1])
	        elif sys.argv[argv_i].split(':')[0] == 'MinLoop' :
	            min_loop = int(sys.argv[argv_i].split(':')[1])
	        elif sys.argv[argv_i].split(':')[0] == 'LoopStability':
	            LS = int(sys.argv[argv_i].split(':')[1])
	        elif sys.argv[argv_i].split(':')[0] == 'StemStability':
	            acceptable_unmatch_num = int(sys.argv[argv_i].split(':')[1])
	        else :
	            print()
	            print('You entered an incorrect option.')
	            print('[MinLength]:[Num]      - Desired minimum length of stem-loop structure')
	            print('[MinComplement]:[Num]  - Desired minimum number of complementary binding pair of stem-loop structure')
	            print('[MinLoop]:[Num]        - Minimum length of loop')
	            print('[LoopStability]:[Num]  - The number of complementary base pairs that can be nested within the loop.')
	            print('[StemStability]:[Num]  - The number of complementary base pairs that can be nested additionally within the stem.')
	            print()
	            ex = 1
	    except :
	        pass
	
	if ex == 1 :
	    sys.exit()
	
	if MSL == 'default' :
	    MSL = 15
	if MCN == 'default' :
	    MCN = MSL - 4
	if min_loop == 'default' :
	    min_loop = 3
	if LS == 'default' :
	    LS = 0
	if acceptable_unmatch_num == 'default' :
	    acceptable_unmatch_num = 0
	
	
	
	SRN = open("./{}_stemloop_repeat_num_7556.txt".format(sys.argv[1][0:-4]),'a')
	SRN.close()
	
	SRN = ("./{}_num_7556.txt".format(sys.argv[1][0:-4]))
	
	OFN = []
	for output_file_num in open("./{}_stemloop_repeat_num_7556.txt".format(sys.argv[1][0:-4]), 'r') :
	    OFN.append(output_file_num)
	
	read_file = ''
	if OFN == [] :
	    read_file = sys.argv[1]
	else :
	    read_file = ("./{}_stemloop_{}.txt".format(sys.argv[1][0:-4], len(OFN)))
	
	
	# def Alignmenr_Match : 두 seq를 받아, alignment를 수행한다.
	"""
	    단순하게, 두 seq가 상보적이면 1 점, 그렇지 않으면 -inf 점을 부여한다. Match or Gap을 유지시킨다. gap score는 -1 점이다. 
	1. matrix인 M과 backtrace를 수행하기 위한 Mdir을 생성한다. 
	2. M과 Mdir의 iniation을 수행한다. - (i,0) & (0,j) 채우기
	3. M을 채운다 : 상보적 쌍일 경우 1 점, 상보적이지 않을 경우 쌍을 이루는 것이 불가능하게 -inf 점을 부여한다. 
	4. Mdir을 채운다 : 이때 M으로 가는 방향의 점수가 같을 경우, random하게 방향이 결정된다. <- 나중에 수정할 필요가 있겠지. 
	5. Mdir에 trace back을 수행하여, alignment 결과를 반환한다. 
	"""
	
	
	def Alignment_Match(seqX, seqY):
	
	    # 1
	    M = []
	    Mdir = []
	    g_score = - 1
	    for i in range(0, len(seqX) + 1):
	        M.append([])
	        Mdir.append([])
	        for j in range(0, len(seqY) + 1):
	            M[i].append([])
	            Mdir[i].append([])
	
	    # 2
	    M[0][0] = 0
	    Mdir[0][0] = 'E'
	    for i in range(1, len(seqX) + 1):
	        M[i][0] = g_score * i
	        Mdir[i][0] = 'u'
	        for j in range(1, len(seqY) + 1):
	            M[0][j] = g_score * j
	            Mdir[0][j] = 'l'
	
	    # 3
	    for i in range(1, len(seqX) + 1):
	        for j in range(1, len(seqY) + 1):
	            if (seqX[i - 1] == 'A' and seqY[j - 1] == 'T') or (seqX[i - 1] == 'T' and seqY[j - 1] == 'A') or (
	                    seqX[i - 1] == 'G' and seqY[j - 1] == 'C') or (seqX[i - 1] == 'C' and seqY[j - 1] == 'G'):
	                m_score = 1
	            elif (seqX[i-1] == '@') or (seqY[j-1] == '@') :
	                m_score = 0 # @는 재실행할 때 사용되는 것. 그래서 없는 취급할 겸 0점
	            else:
	                m_score = float('inf') * (-1)
	            M[i][j] = max(M[i - 1][j - 1] + m_score, M[i - 1][j] + g_score, M[i][j - 1] + g_score)
	
	            mdir = ''
	            if M[i][j] == M[i - 1][j - 1] + m_score:
	                mdir += 'd'
	            if M[i][j] == M[i - 1][j] + g_score:
	                mdir += 'u'
	            if M[i][j] == M[i][j - 1] + g_score:
	                mdir += 'l'
	            Mdir[i][j] = mdir
	
	    # 4
	    sequenceX = ''
	    sequenceY = ''
	
	    num_match = 0
	    while i != 0 or j != 0:
	        if 'd' in Mdir[i][j]:
	            cur_dir = 'd'
	        elif 'u' in Mdir[i][j]:
	            cur_dir = 'u'
	        elif 'l' in Mdir[i][j]:
	            cur_dir = 'l'
	
	        if cur_dir == 'd':
	            # cur_dir = random.choice(Mdir[i-1][j-1])
	            sequenceX += seqX[i - 1]
	            sequenceY += seqY[j - 1]
	            i = i - 1
	            j = j - 1
	            num_match += 1
	        elif cur_dir == 'u':
	            # cur_dir = random.choice(Mdir[i-1][j])
	            sequenceX += seqX[i - 1]
	            sequenceY += '-'
	            i = i - 1
	        elif cur_dir == 'l':
	            # cur_dir = random.choice(Mdir[i][j-1])
	            sequenceX += '-'
	            sequenceY += seqY[j - 1]
	            j = j - 1
	
	    # 5
	    return_list = []
	    return_list.append(num_match)
	    return_list.append(sequenceX[::-1])
	    return_list.append(sequenceY[::-1])
	
	    return return_list
	
	
	# stem_lopp를 가질 가능성이 높은 파트를 찾는 탐색 단계
	"""
	1. seq_name을 저장한다.
	2. threefold_line이 항상 3줄(line, f.line, ff.line)을 가지도록 만든다. 
	    2-1. 초반에는 threefold_line이 3줄 모두 채워져 있지 않다. 때문에 우선 threefold_line을 3줄 채운다. 
	    2-2. 초반 이후부터는 : 가장 오래된 line을 지우고, 새로운 line을 받는다. 이 과정으로 threefold_line을 계속 유지시킨다. 
	3. threefold_line에서 step이 움직인다. step의 앞 seq와 뒷 seq로 alginment를 진행한다. 
	4. alignment 결과가 MCN 기준을 만족할 경우, 해당 locate를 pre_locate_list에 저장한다. 
	5. 마지막에는 겸사겸사 seq_length도 구한다. 
	결과 : 우리는 이제 stem-loop구조를 가질 가능성이 높은 ( = MCN 기준을 만족한) 부위들을, 파일 내 seq로부터 추출하였다. 
	    해당 정보는 pre_locate_list에 저장되어 있다. 
	"""
	
	"""
	기본 원리 : 
	너무 긴 seq를 받을시, 메모리에 부담을 주지 않기 위해, 
	파일을 현재 한 줄, 앞의 한 줄, 뒤의 한 줄씩, 총 3줄을 가지고 stem-loop를 찾는다. 
	"""
	# MSL = int(sys.argv[2])       # MSL, argv[2] : 최소의 stem-loop 구조 길이. (Minimum Structure Length)
	# MCN = int(sys.argv[3])       # MCN, argv[3] : MSL에서 제작된 stem이 가져야 하는 최소의 상보쌍 결합 수 (Minimum Complementary pair Number)
	
	
	seq_name = ''                # seq_name : input으로 받은 파일의 이름을 담기위한 변수.
	
	line_num = 1                 # line_num : 현재 읽고 있는 line의 파일 내 위치
	threefold_line = []          # threefold_line : 현재 line, 이전 line, 이후 line, 총 3줄을 담기위한 리스트
	lines = ''                   # lines : threefold_line (list)의 각 요소들을 하나로 합친 str
	
	step = MSL + 0               # step : alignment를 수행하기위한 lines에서의 현재 위치.
	"""
	lines에서 alignment를 수행할 때, step을 기준으로 앞의 seq와 뒤의 seq를 가지고 alignment를 수행한다. 
	이때 step이 file의 1번째부터 위치하게 되면, step 앞의 seq는 만들어지지 못해, range error가 발생한다. 
	때문에 step의 초기값은 MSL + 0 이다. 
	"""
	
	locate = copy.deepcopy(MSL)  # locate : 파일의 seq 내에서 현재 위치. step이 MSL+0에서부터 시작하므로, locate 또한 MSL+0에서부터 시작한다.
	pre_locate_list = []         # pre_locate_list : MCN 기준을 충족한 stem-loop구조는, 그 locate를 pre_locate_list에 저장해둔다.
	seq_length = ''              # seq_length : 파일 내 seq의 총 길이
	
	# line, f.line, ff.line : threefold_line을 이루는 각 line들. 중요 : < step은 f.line 내에서 움직인다. >
	line = ''
	front_line = ''
	front_front_line = ''
	
	
	for line in open(read_file):
	    line = line.rstrip()
	    if line[0] == '>':
	        seq_name = line
	    elif line.startswith('#') :
	        break
	    else:
	        if line_num <= 3:
	            threefold_line.append(line)
	
	            # threefold - Begining part
	            if line_num == 3:
	                lines = ''.join(threefold_line)
	                while step <= len(front_front_line) + len(front_line) - 1:
	                    Left = lines[step - MSL: step]
	                    Right = lines[step + MSL - 1: step - 1: -1]
	                    if MCN <= Alignment_Match(Left, Right)[0]:
	                        pre_locate_list.append(locate)
	
	                    step += 1
	                    locate += 1
	
	            line_num += 1
	
	        # threefold - Middle part
	        else:
	            del (threefold_line[0])
	            threefold_line.append(line)
	            lines = ''.join(threefold_line)
	
	            step = len(front_line)
	            while step <= len(front_front_line) + len(front_line) - 1:
	                Left = lines[step - MSL: step]
	                Right = lines[step + MSL - 1: step - 1: -1]
	                if MCN <= Alignment_Match(Left, Right)[0]:
	                    pre_locate_list.append(locate)
	
	                step += 1
	                locate += 1
	            step = len(front_line)
	
	        front_front_line = copy.deepcopy(front_line)
	        front_line = copy.deepcopy(line)
	
	# threefold - Last part
	step = len(lines) - len(front_line)
	while step <= len(lines) - 1:
	    Left = lines[step - MSL: step]
	    Right = lines[step + MSL - 1: step - 1: -1]
	    if MCN <= Alignment_Match(Left, Right)[0]:
	        pre_locate_list.append(locate)
	
	    step += 1
	    locate += 1
	
	seq_length = copy.deepcopy(locate)
	
	
	# pre_locate_list 값 정리.
	"""
	pre_locate_list에서는, 한 칸 한 칸 앞으로 전진하는 step 특성상, 하나의 stem-loop구조를 여러 요소들이 지정하는 경우가 빈번하다. 
	가령 MSL을 15로 지정할 경우, 길이가 20인 stem-loop구조를 찾게된다면, 해당 stem-loop를 찾은 step들은 약 5개 정도일 것이다. 
	때문에, pre_locate_list 내에서는 하나의 stem-lopp를 가리키는 중복된 요소들이 다수 존재하며, 이 요소들은 연속되어 있는 경우가 대다수이다. 
	이에 연속된 값들 중 하나만을 남겨놓고, 모두 제거하는 작업을 거친다. 
	"""
	for lo_num in range(1, len(pre_locate_list)):
	    if pre_locate_list[lo_num] == pre_locate_list[lo_num - 1] + 1:
	        pre_locate_list[lo_num - 1] = 'X'
	while 'X' in pre_locate_list:
	    pre_locate_list.remove('X')
	
	
	
	# 형성 단계 : 탐색된 결과를 토대로, stem-loop구조를 형성시킨다.
	"""
	형성 단계는 탐색 단계와 달리, stem-loop를 확장시킨다. 
	loop구조를 형성화하며, 
	보다 실제 stem-loop 구조에 가깝도록, stem을 늘이거나 줄여, stem의 길이를 보정한다. 
	argv[4]가 최소한의 loop 길이를 설정하며, 
	argv[5]가 loop 내에 존재할 수 있는 상보적 염기쌍의 수를 정의한다. 
	"""
	# min_loop = int(sys.argv[4])  # min_loop, argv[4] : 최소 loop의 길이
	
	
	"""
	초기 : 1 line과 2 line을 합친 lines에서, step이 1 line에서 활동한다. 
	중기 : (1,2,3) lines ~ (n-2,n-1,n) lines에서, step이 각각 2 ~ n-1 line에서 활동한다. 
	후기 : n-1,n lines에서, step이 n line에서 활동한다. 
	"""
	
	loop = ''                    # loop : loop seq를 담는 변수
	
	loop_locate_map = []         # loop_locate_map :
	
	# line, f.line, ff.line : threefold_line을 이루는 각 line들.
	line = ''
	front_line = ''
	front_front_line = ''
	
	line_num = 1                 # line_num : 파일 내 seq의 line 위치
	threefold_line = []          # threefold_line : line, f.line, ff.line을 담는 list
	lines = ''                   # lines : threefold_line을 합친 str
	
	step = ''                    # step : threefold_line 내에서의 위치. step을 기준으로 loop가 형성되며, alignment가 시행된다.
	locate = ''                  # locate : 파일 내 seq에서의 위치.
	
	sum_len_line = 0             # sum_len_line :
	ll_num = 0                   # ll_num :
	
	
	# loop_stability : loop 구조의 형성 조절.
	"""
	loop의 크기는 환경에 영향을 받는다. 
	간단히는, loop의 크기는 비상보적 염기쌍만큼의 크기를 갖겠지만, 
	상보적 염기쌍이 끊어지기 쉬운 환경에서는, loop 주변의 stem에서 염기쌍이 끊어져, loop의 크기가 더 커질 수 있을 것이다. 
	이때, loop_stability는 loop 주변의 stem에서, 상보적 염기쌍이 끊어질 수 있는 갯수를 설정한 값이다. 
	만약 TTTTCAAACAAAA seq에서 
	loop_stability가 0이라면, loop는 CAAAC일 것이며, 
	loop_stability가 1이라면, loop는 TCAAACA일 것이다.  
	"""
	# loop_stability = int(sys.argv[5])
	
	for line in open(read_file):
	    line = line.rstrip()
	    if line.startswith('#') :
	        break
	    if line[0] != '>' :
	        # second rotate - begining part
	        if line_num == 3:
	            while 1:
	                try :
	                    if pre_locate_list[ll_num] <= len(front_front_line) - 1:
	                        locate = pre_locate_list[ll_num]
	                        step = pre_locate_list[ll_num]
	
	                        Left = lines[step - MSL: step]
	                        Right = lines[step + MSL - 1: step - 1: -1]
	                        AM_list = Alignment_Match(Left, Right)
	
	                        loop_stability = LS
	
	                        re_num = 0
	                        loop = ''
	                        left_loop_locate = 0
	                        right_loop_locate = 0
	
	                        while (loop_stability != -1) and (re_num != len(AM_list[1])):
	                            if AM_list[1][len(AM_list[1]) - 1 - re_num] == '-':
	                                loop = loop + AM_list[2][len(AM_list[1]) - 1 - re_num]
	                                right_loop_locate += 1
	                            elif (AM_list[2][len(AM_list[1]) - 1 - re_num] == '-'):
	                                loop = AM_list[1][len(AM_list[1]) - 1 - re_num] + loop
	                                left_loop_locate += 1
	                            else:
	                                if loop_stability != 0:
	                                    loop = AM_list[1][len(AM_list[1]) - 1 - re_num] + loop + AM_list[2][
	                                        len(AM_list[1]) - 1 - re_num]
	                                    left_loop_locate += 1
	                                    right_loop_locate += 1
	                                loop_stability -= 1
	                            re_num += 1
	
	                        if len(loop) >= min_loop:
	                            loop_locate_list = []
	                            if loop == '':
	                                loop = 'no loop'
	                            loop_locate_list.append(loop)
	                            loop_locate_list.append('loop locate : ')
	                            loop_locate_list.append(locate - left_loop_locate)
	                            loop_locate_list.append(locate + right_loop_locate - 1)
	                            loop_locate_map.append(loop_locate_list)
	
	
	                        ll_num += 1
	
	                    else:
	                        break
	
	                except IndexError:
	                    break
	
	        if line_num <= 3:
	            threefold_line.append(line)
	            line_num += 1
	            lines = ''.join(threefold_line)
	
	        if line_num > 3:
	            if line_num == 5:
	                del (threefold_line[0])
	                threefold_line.append(line)
	                lines = ''.join(threefold_line)
	            else:
	                line_num += 1
	
	            # second rotate - middle part
	            while 1:
	                try:
	                    if sum_len_line - 1 < pre_locate_list[ll_num] <= sum_len_line + len(front_line) - 1:
	                        locate = pre_locate_list[ll_num]
	                        step = len(front_front_line) + (locate - sum_len_line)
	
	                        Left = lines[step - MSL: step]
	                        Right = lines[step + MSL - 1: step - 1: -1]
	                        AM_list = Alignment_Match(Left, Right)
	                        loop_stability = LS
	
	                        re_num = 0
	                        loop = ''
	                        left_loop_locate = 0
	                        right_loop_locate = 0
	
	                        while (loop_stability != -1) and (re_num != len(AM_list[1])):
	                            if AM_list[1][len(AM_list[1]) - 1 - re_num] == '-':
	                                loop = loop + AM_list[2][len(AM_list[1]) - 1 - re_num]
	                                right_loop_locate += 1
	                            elif (AM_list[2][len(AM_list[1]) - 1 - re_num] == '-'):
	                                loop = AM_list[1][len(AM_list[1]) - 1 - re_num] + loop
	                                left_loop_locate += 1
	                            else:
	                                if loop_stability != 0:
	                                    loop = AM_list[1][len(AM_list[1]) - 1 - re_num] + loop + AM_list[2][
	                                        len(AM_list[1]) - 1 - re_num]
	                                    left_loop_locate += 1
	                                    right_loop_locate += 1
	                                loop_stability -= 1
	                            re_num += 1
	
	                        if len(loop) >= min_loop:
	                            loop_locate_list = []
	                            if loop == '' :
	                                loop = 'no loop'
	                            loop_locate_list.append(loop)
	                            loop_locate_list.append('loop locate : ')
	                            loop_locate_list.append(locate - left_loop_locate)
	                            loop_locate_list.append(locate + right_loop_locate - 1)
	                            loop_locate_map.append(loop_locate_list)
	
	                        ll_num += 1
	
	                    else:
	                        break
	
	                except IndexError:
	                    break
	
	        front_front_line = copy.deepcopy(front_line)
	        front_line = copy.deepcopy(line)
	
	        sum_len_line += len(front_front_line)
	
	# second rotate - last part
	while 1:
	    try:
	        if sum_len_line + len(front_line) - 1 < pre_locate_list[ll_num]:
	            locate = pre_locate_list[ll_num]
	            step = len(front_front_line) + (locate - sum_len_line)
	
	            Left = lines[step - MSL: step]
	            Right = lines[step + MSL - 1: step - 1: -1]
	            AM_list = Alignment_Match(Left, Right)
	            loop_stability = LS
	
	            re_num = 0
	            loop = ''
	            left_loop_locate = 0
	            right_loop_locate = 0
	
	            while (loop_stability != -1) and (re_num != len(AM_list[1])):
	                if AM_list[1][len(AM_list[1]) - 1 - re_num] == '-':
	                    loop = loop + AM_list[2][len(AM_list[1]) - 1 - re_num]
	                    right_loop_locate += 1
	                elif (AM_list[2][len(AM_list[1]) - 1 - re_num] == '-'):
	                    loop = AM_list[1][len(AM_list[1]) - 1 - re_num] + loop
	                    left_loop_locate += 1
	                else:
	                    if loop_stability != 0:
	                        loop = AM_list[1][len(AM_list[1]) - 1 - re_num] + loop + AM_list[2][
	                            len(AM_list[1]) - 1 - re_num]
	                        left_loop_locate += 1
	                        right_loop_locate += 1
	                    loop_stability -= 1
	                re_num += 1
	
	            if len(loop) >= min_loop:
	                loop_locate_list = []
	                if loop == '':
	                    loop = 'no loop'
	                loop_locate_list.append(loop)
	                loop_locate_list.append('loop locate : ')
	                loop_locate_list.append(locate - left_loop_locate)
	                loop_locate_list.append(locate + right_loop_locate - 1)
	                loop_locate_map.append(loop_locate_list)
	
	            ll_num += 1
	
	        else:
	            break
	
	    except IndexError:
	        break
	
	
	# stem 구체화 단계 : 제작된 stem-loop구조를 토대로, stem의 길이를 조정한다.
	"""
	앞의 loop 보정과 마찬가지로, stem의 길이가 더 늘어날 수 있는가에 대해 진행되는 조정 작업이다. 
	환경의 변화로, stem이 비상보적인 염기쌍을 더 내포할 수 있다면, stem의 길이는 더 길어질 것이다. 
	stem-loop 구조의 5'end 와 3'end 를 더 확장시켜, 
	acceptable_unmatch_num (argv[6])에서 지정된 값만큼, stem 내에서의 비상보적 염기쌍을 추가로 허용해준다. 
	"""
	# acceptable_unmatch_num = int(sys.argv[6])
	
	## third rotate
	correction = 0
	
	stem_map = []
	
	line_num = 1
	threefold_line = []
	
	line = ''
	front_line = ''
	front_front_line = ''
	lines = ''
	
	sum_len_line = 0
	ll_num = 0
	
	for line in open(read_file) :
	    line = line.rstrip()
	    if line.startswith('#') :
	        break
	    if line[0] != '>' :
	        # third rotate - begining part
	        if line_num == 3:
	            while 1:
	                if ll_num < len(loop_locate_map):
	                    if loop_locate_map[ll_num][3] <= len(front_front_line) - 1:
	
	                        correction = 0
	                        try:
	                            while 1:
	                                Left = lines[loop_locate_map[ll_num][2] - (
	                                        2 * MSL - len(loop_locate_map[ll_num][0])) // 2 - correction :
	                                             loop_locate_map[ll_num][2]]
	                                Right = lines[loop_locate_map[ll_num][3] + (
	                                        2 * MSL - len(loop_locate_map[ll_num][0])) // 2 + correction :
	                                              loop_locate_map[ll_num][3]: -1]
	                                AM_list = Alignment_Match(Left, Right)
	
	                                s_num = 0
	                                while (AM_list[1][s_num] == '-') or (AM_list[2][s_num] == '-'):
	                                    s_num += 1
	
	                                if s_num > acceptable_unmatch_num:
	                                    stem_L = AM_list[1][s_num:]
	                                    stem_R = AM_list[2][s_num:]
	                                    stem_R = stem_R[::-1]
	                                    stem_list = [stem_L, stem_R]
	                                    stem_map.append(stem_list)
	                                    break
	                                else:
	                                    correction += 1
	                        except:
	                            stem_L = AM_list[1][s_num:]
	                            stem_R = AM_list[2][s_num:]
	                            stem_R = stem_R[::-1]
	                            stem_list = [stem_L, stem_R]
	                            stem_map.append(stem_list)
	
	                        ll_num += 1
	
	                    else:
	                        break
	
	                else:
	                    break
	
	        if line_num <= 3:
	            threefold_line.append(line)
	            line_num += 1
	            lines = ''.join(threefold_line)
	
	        if line_num > 3:
	            if line_num == 5:
	                del (threefold_line[0])
	                threefold_line.append(line)
	                lines = ''.join(threefold_line)
	            else:
	                line_num += 1
	
	            # third rotate - middle part
	            while 1:
	                if ll_num < len(loop_locate_map):
	                    if sum_len_line - 1 < loop_locate_map[ll_num][3] <= sum_len_line + len(front_line) - 1:
	
	                        correction = 0
	                        try:
	                            while 1:
	                                Left = lines[loop_locate_map[ll_num][2] - (
	                                        2 * MSL - len(loop_locate_map[ll_num][0])) // 2 - correction - (
	                                                     sum_len_line - len(front_front_line)): loop_locate_map[ll_num][
	                                                                                                2] - (
	                                                                                                    sum_len_line - len(
	                                                                                                front_front_line))]
	                                Right = lines[loop_locate_map[ll_num][3] + (
	                                        2 * MSL - len(loop_locate_map[ll_num][0])) // 2 + correction - (
	                                                      sum_len_line - len(front_front_line)):
	                                              loop_locate_map[ll_num][3] - (sum_len_line - len(front_front_line)): -1]
	                                AM_list = Alignment_Match(Left, Right)
	
	                                s_num = 0
	                                while (AM_list[1][s_num] == '-') or (AM_list[2][s_num] == '-'):
	                                    s_num += 1
	
	                                if s_num > acceptable_unmatch_num:
	                                    stem_L = AM_list[1][s_num:]
	                                    stem_R = AM_list[2][s_num:]
	                                    stem_R = stem_R[::-1]
	                                    stem_list = [stem_L, stem_R]
	                                    stem_map.append(stem_list)
	                                    break
	                                else:
	                                    correction += 1
	                        except:
	                            stem_L = AM_list[1][s_num:]
	                            stem_R = AM_list[2][s_num:]
	                            stem_R = stem_R[::-1]
	                            stem_list = [stem_L, stem_R]
	                            stem_map.append(stem_list)
	
	                        ll_num += 1
	
	                    else:
	                        break
	
	                else:
	                    break
	
	        front_front_line = copy.deepcopy(front_line)
	        front_line = copy.deepcopy(line)
	
	        sum_len_line += len(front_front_line)
	
	# third rotate - last part
	while 1:
	    if ll_num < len(loop_locate_map):
	        if sum_len_line + len(front_line) - 1 < loop_locate_map[ll_num][3]:
	
	            correction = 0
	            try:
	                while 1:
	                    Left = lines[loop_locate_map[ll_num][2] - (
	                            2 * MSL - len(loop_locate_map[ll_num][0])) // 2 - correction - (
	                                         sum_len_line - len(front_front_line)): loop_locate_map[ll_num][2] - (
	                            sum_len_line - len(front_front_line))]
	                    Right = lines[loop_locate_map[ll_num][3] + (
	                            2 * MSL - len(loop_locate_map[ll_num][0])) // 2 + correction - (
	                                          sum_len_line - len(front_front_line)): loop_locate_map[ll_num][3] - (
	                            sum_len_line - len(front_front_line)): -1]
	                    AM_list = Alignment_Match(Left, Right)
	
	                    s_num = 0
	                    while (AM_list[1][s_num] == '-') or (AM_list[2][s_num] == '-'):
	                        s_num += 1
	
	                    if s_num > acceptable_unmatch_num:
	                        stem_L = AM_list[1][s_num:]
	                        stem_R = AM_list[2][s_num:]
	                        stem_R = stem_R[::-1]
	                        stem_list = [stem_L, stem_R]
	                        stem_map.append(stem_list)
	                        break
	                    else:
	                        correction += 1
	            except:
	                stem_L = AM_list[1][s_num:]
	                stem_R = AM_list[2][s_num:]
	                stem_R = stem_R[::-1]
	                stem_list = [stem_L, stem_R]
	                stem_map.append(stem_list)
	
	            ll_num += 1
	
	        else:
	            break
	
	    else:
	        break
	
	
	
	num = 1
	true_final = []
	for aaaa in range(0, len(loop_locate_map)):
	    num += 1
	    sa_num = 0
	    sb_num = 0
	    for sa in stem_map[aaaa][0] :
	        if sa == '-' :
	            sa_num += 1
	    for sb in stem_map[aaaa][1] :
	        if sb == '-' :
	            sb_num += 1
	    true_final.append([loop_locate_map[aaaa][2] + 1 - (len(stem_map[aaaa][0])-sa_num),
	        loop_locate_map[aaaa][2] + 1, loop_locate_map[aaaa][3] + 1,
	        loop_locate_map[aaaa][3] + 1 + (len(stem_map[aaaa][1])-sb_num)])
	
	
	
	
	
	# 재실행
	true_final.append([-1,-1,-1,-1])
	output = open("./output.txt",'w')
	locate = 1                   # locate : 파일의 seq 내에서 현재 위치.
	i = 0
	letter_num = 0
	for line in open(read_file):
	    line = line.rstrip()
	    if line.startswith('#') :
	        break
	    if line[0] != '>' :
	        for letter in line :
	            if true_final[i][0] <= locate <= true_final[i][3] :
	                if locate == true_final[i][3]:
	                    output.write('@')
	                    letter_num += 1
	                    i += 1
	            else :
	                output.write(letter)
	                letter_num += 1
	
	            locate += 1
	
	            if letter_num == 70 :
	                output.write('\n')
	                letter_num = 0
	
	
	if not true_final == [] :
	    output.write('\n')
	    for i in true_final :
	        output.write('#')
	        output.write(str(i[0]))
	        output.write('/')
	        output.write(str(i[1]))
	        output.write('/')
	        output.write(str(i[2]))
	        output.write('/')
	        output.write(str(i[3]))
	        output.write('\n')
	
	output.close()
	
	
	SRN = open("./{}_stemloop_repeat_num_7556.txt".format(sys.argv[1][0:-4]),'a')
	for i in range(0,len(true_final)-1) :
	    for j in range(0,4) :
	        SRN.write(str(true_final[i][j]))
	        SRN.write('/')
	    SRN.write(';')
	SRN.write('\n')
	SRN.close()
	
	OFN = []
	for output_file_num in open("./{}_stemloop_repeat_num_7556.txt".format(sys.argv[1][0:-4]), 'r') :
	    OFN.append(output_file_num)
	
	os.rename("./output.txt", "./{}_stemloop_{}.txt".format(sys.argv[1][0:-4], len(OFN)))
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	print()
	print(
	    '=============================================================================================================================================')
	print(' NAME :', seq_name)
	print(' simple range : ', '(', 1, '...', seq_length, ')')
	print()
	print(
	    '=============================================================================================================================================')
	print(' num\t', 'loop\t\t', 'stem')
	print(
	    '=============================================================================================================================================')
	num = 1
	true_final = []
	for aaaa in range(0, len(loop_locate_map)):
	    print(' ', num, end=')\t')
	    num += 1
	    print(loop_locate_map[aaaa][0])
	    sa_num = 0
	    sb_num = 0
	
	    for sa in stem_map[aaaa][0] :
	        if sa == ('-' or '@') :
	            sa_num += 1
	    for sb in stem_map[aaaa][1] :
	        if sb == ('-' or '@') :
	            sb_num += 1
	    true_final.append([loop_locate_map[aaaa][2] + 1 - (len(stem_map[aaaa][0])-sa_num),
	        loop_locate_map[aaaa][2] + 1, loop_locate_map[aaaa][3] + 1,
	        loop_locate_map[aaaa][3] + 1 + (len(stem_map[aaaa][1])-sb_num)])
	    print("     5'", stem_map[aaaa][0], "3'")
	    print("     3'", stem_map[aaaa][1][::-1], "5'")
	    print(
	        '---------------------------------------------------------------------------------------------------------------------------------------------')
	
	print()
	print()
	
	
	
	
	
	
	
	if not loop_locate_map == [] :
	    os.execl(sys.executable, sys.executable, *sys.argv)
	else :
	    os.remove("./{}_stemloop_{}.txt".format(sys.argv[1][0:-4], len(OFN)))
	
	    zero_list = []
	    zero_list_list = []
	    for i in range(1, len(OFN)) :
	        for line in open("./{}_stemloop_{}.txt".format(sys.argv[1][0:-4], i)) :
	            if line.startswith('#') :
	                line = line.rstrip()
	                line = line[1:]
	                line = line.split('/')
	                zero_list_list.append(line)
	        zero_list.append(zero_list_list)
	        zero_list_list = []
	
	    zero = open(("./{}_stemloop_result_0.txt".format(sys.argv[1][0:-4])), 'a')
	    zero.close()
	
	    file1 = ('./{}'.format(sys.argv[1]))
	    file2 = ("./{}_stemloop_result_0.txt".format(sys.argv[1][0:-4]))
	    shutil.copy(file1, file2)
	
	
	
	    letter_num = 0
	    h = 0
	    i = 0
	    j = 0
	    k = 0
	    while h != len(zero_list) :
	        num70 = 0
	        pass_switch = 0
	        result = open("./{}_stemloop_result_{}.txt".format(sys.argv[1][0:-4], k+1), 'a')
	        for line in open("./{}_stemloop_result_{}.txt".format(sys.argv[1][0:-4], h), 'r') :
	            line = line.rstrip()
	
	            if line.startswith('>') :
	                result.write(line)
	                result.write('\n')
	
	            else :
	                for letter in line :
	                    num70 += 1
	                    if h == len(zero_list):
	                        pass
	                    else :
	                        if (letter == 'B') or (letter == 'O') or (letter == 'H') or (letter == 'D') :
	                            pass_switch = 1
	                            letter_num += 1
	
	
	                        if pass_switch == 0 :
	                            letter_num += 1
	                            if letter_num == int(zero_list[h][i][j]) :
	                                if (letter == 'A') and (j == 0) :
	                                    letter = 'B'
	                                elif (letter == 'A') and (j == 3) :
	                                    letter = 'b'
	                                elif (letter == 'A') and (j == 1):
	                                    letter = 'E'
	                                elif (letter == 'A') and (j == 2):
	                                    letter = 'e'
	                                elif (letter == 'T') and (j == 0) :
	                                    letter = 'O'
	                                elif (letter == 'T') and (j == 3) :
	                                    letter = 'o'
	                                elif (letter == 'T') and (j == 1):
	                                    letter = 'P'
	                                elif (letter == 'T') and (j == 2):
	                                    letter = 'p'
	                                elif (letter == 'G') and (j == 0) :
	                                    letter = 'H'
	                                elif (letter == 'G') and (j == 3) :
	                                    letter = 'h'
	                                elif (letter == 'G') and (j == 1):
	                                    letter = 'I'
	                                elif (letter == 'G') and (j == 2):
	                                    letter = 'i'
	                                elif (letter == 'C') and (j == 0) :
	                                    letter = 'D'
	                                elif (letter == 'C') and (j == 3) :
	                                    letter = 'd'
	                                elif (letter == 'C') and (j == 1):
	                                    letter = 'F'
	                                elif (letter == 'C') and (j == 2):
	                                    letter = 'f'
	
	                                j += 1
	                                if j == 4 :
	                                    i += 1
	                                    j = 0
	                                    if zero_list[h][i][j] == '-1' :
	                                        h += 1
	                                        i = 0
	
	                        if (letter == 'b') or (letter == 'o') or (letter == 'h') or (letter == 'd') :
	                            pass_switch = 0
	
	                    result.write(letter)
	
	                    if num70 % 70 == 0 :
	                        result.write('\n')
	        letter_num = 0
	        k += 1
	        result.close()
	
	
	
	    true_true_true_final = open("./{}_stemloop_final_result.txt".format(sys.argv[1][0:-4]), 'a')
	    for line in open("./{}_stemloop_result_{}.txt".format(sys.argv[1][0:-4], k), 'r') :
	        if line.startswith('>') :
	            true_true_true_final.write(line)
	        else :
	            line = line.replace('B','(A')
	            line = line.replace('b','A)')
	            line = line.replace('E','[A')
	            line = line.replace('e','A]')
	            line = line.replace('O','(T')
	            line = line.replace('o','T)')
	            line = line.replace('P','[T')
	            line = line.replace('p','T]')
	            line = line.replace('H','(G')
	            line = line.replace('h','G)')
	            line = line.replace('I','[G')
	            line = line.replace('i','G]')
	            line = line.replace('D','(C')
	            line = line.replace('d','C)')
	            line = line.replace('F','[C')
	            line = line.replace('f','C]')
	
	            true_true_true_final.write(line)
	    true_true_true_final.close()
	

if __name__ == '__main__':
	RNAmap()

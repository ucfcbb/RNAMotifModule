import math

def get_initialized_basepair_index():

	n = 289
	basepair_index = [0 for i in range(n)]

	#	<*****matrix_0*****>
	basepair_index[0] = 0

	#	the theoretically impossible base pairs are treated as hydrogen bonds
	#	Note that if basepair_index[i] > 0, then it means it is a possible base pair
	#	the score should be located at (basepair_index[i] - 1) of the isosteric base pair
	#	substitution matrix


	#	<*****matrix_1*****>
	basepair_index[1] = 4		#	cis	W/W	AA	(0)
	basepair_index[2] = 2		#	cis	W/W	AC	(1)
	basepair_index[3] = 3		#	cis	W/W	AG	(2)
	basepair_index[4] = 1		#	cis	W/W	AU	(3)
	basepair_index[5] = 2		#	cis	W/W	CA	(4)
	basepair_index[6] = 6		#	cis	W/W	CC	(5)
	basepair_index[7] = 1		#	cis	W/W	CG	(6)
	basepair_index[8] = 5		#	cis	W/W	CU	(7)
	basepair_index[9] = 3		#	cis	W/W	GA	(8)
	basepair_index[10] = 1		#	cis	W/W	GC	(9)
	basepair_index[11] = 0		#	cis	W/W	GG	(10)
	basepair_index[12] = 2		#	cis	W/W	GU	(11)
	basepair_index[13] = 1		#	cis	W/W	UA	(12)
	basepair_index[14] = 5		#	cis	W/W	UC	(13)
	basepair_index[15] = 2		#	cis	W/W	UG	(14)
	basepair_index[16] = 6		#	cis	W/W	UU	(15)

	#	<*****matrix_2*****>
	basepair_index[17] = 10	#	trans	W/W	AA	(16)
	basepair_index[18] = 9		#	trans	W/W	AC	(17)
	basepair_index[19] = 0		#	trans	W/W	AG	(18)
	basepair_index[20] = 7		#	trans	W/W	AU	(19)
	basepair_index[21] = 9		#	trans	W/W	CA	(20)
	basepair_index[22] = 12	#	trans	W/W	CC	(21)
	basepair_index[23] = 8		#	trans	W/W	CG	(22)
	basepair_index[24] = 11	#	trans	W/W	CU	(23)
	basepair_index[25] = 0		#	trans	W/W	GA	(24)
	basepair_index[26] = 8		#	trans	W/W	GC	(25)
	basepair_index[27] = 10	#	trans	W/W	GG	(26)
	basepair_index[28] = 9		#	trans	W/W	GU	(27)
	basepair_index[29] = 7		#	trans	W/W	UA	(28)
	basepair_index[30] = 11	#	trans	W/W	UC	(29)
	basepair_index[31] = 9		#	trans	W/W	UG	(30)
	basepair_index[32] = 12	#	trans	W/W	UU	(31)

	#	<*****matrix_3*****>
	basepair_index[33] = 0		#	cis	W/H	AA	(32)
	basepair_index[34] = 0		#	cis	W/H	AC	(33)
	basepair_index[35] = 15	#	cis	W/H	AG	(34)
	basepair_index[36] = 15	#	cis	W/H	AU	(35)
	basepair_index[37] = 0		#	cis	W/H	CA	(36)
	basepair_index[38] = 14	#	cis	W/H	CC	(37)
	basepair_index[39] = 13	#	cis	W/H	CG	(38)
	basepair_index[40] = 13	#	cis	W/H	CU	(39)
	basepair_index[41] = 15	#	cis	W/H	GA	(40)
	basepair_index[42] = 0		#	cis	W/H	GC	(41)
	basepair_index[43] = 16	#	cis	W/H	GG	(42)
	basepair_index[44] = 0		#	cis	W/H	GU	(43)
	basepair_index[45] = 13	#	cis	W/H	UA	(44)
	basepair_index[46] = 0		#	cis	W/H	UC	(45)
	basepair_index[47] = 13	#	cis	W/H	UG	(46)
	basepair_index[48] = 14	#	cis	W/H	UU	(47)

	#	<*****matrix_4*****>
	basepair_index[49] = 20	#	trans	W/H	AA	(48)
	basepair_index[50] = 0		#	trans	W/H	AC	(49)
	basepair_index[51] = 20	#	trans	W/H	AG	(50)
	basepair_index[52] = 0		#	trans	W/H	AU	(51)
	basepair_index[53] = 18	#	trans	W/H	CA	(52)
	basepair_index[54] = 17	#	trans	W/H	CC	(53)
	basepair_index[55] = 18	#	trans	W/H	CG	(54)
	basepair_index[56] = 0		#	trans	W/H	CU	(55)
	basepair_index[57] = 0		#	trans	W/H	GA	(56)
	basepair_index[58] = 0		#	trans	W/H	GC	(57)
	basepair_index[59] = 21	#	trans	W/H	GG	(58)
	basepair_index[60] = 20	#	trans	W/H	GU	(59)
	basepair_index[61] = 17	#	trans	W/H	UA	(60)
	basepair_index[62] = 0		#	trans	W/H	UC	(61)
	basepair_index[63] = 19	#	trans	W/H	UG	(62)
	basepair_index[64] = 18	#	trans	W/H	UU	(63)

	#	<*****matrix_5*****>
	basepair_index[65] = 22	#	cis	W/S	AA	(64)
	basepair_index[66] = 22	#	cis	W/S	AC	(65)
	basepair_index[67] = 22	#	cis	W/S	AG	(66)
	basepair_index[68] = 22	#	cis	W/S	AU	(67)
	basepair_index[69] = 23	#	cis	W/S	CA	(68)
	basepair_index[70] = 23	#	cis	W/S	CC	(69)
	basepair_index[71] = 23	#	cis	W/S	CG	(70)
	basepair_index[72] = 23	#	cis	W/S	CU	(71)
	basepair_index[73] = 24	#	cis	W/S	GA	(72)
	basepair_index[74] = 24	#	cis	W/S	GC	(73)
	basepair_index[75] = 26	#	cis	W/S	GG	(74)
	basepair_index[76] = 24	#	cis	W/S	GU	(75)
	basepair_index[77] = 25	#	cis	W/S	UA	(76)
	basepair_index[78] = 25	#	cis	W/S	UC	(77)
	basepair_index[79] = 25	#	cis	W/S	UG	(78)
	basepair_index[80] = 25	#	cis	W/S	UU	(79)

	#	<*****matrix_6*****>
	basepair_index[81] = 27	#	trans	W/S	AA	(80)
	basepair_index[82] = 27	#	trans	W/S	AC	(81)
	basepair_index[83] = 27	#	trans	W/S	AG	(82)
	basepair_index[84] = 27	#	trans	W/S	AU	(83)
	basepair_index[85] = 27	#	trans	W/S	CA	(84)
	basepair_index[86] = 27	#	trans	W/S	CC	(85)
	basepair_index[87] = 27	#	trans	W/S	CG	(86)
	basepair_index[88] = 27	#	trans	W/S	CU	(87)
	basepair_index[89] = 0		#	trans	W/S	GA	(88)
	basepair_index[90] = 28	#	trans	W/S	GC	(89)
	basepair_index[91] = 0		#	trans	W/S	GG	(90)
	basepair_index[92] = 28	#	trans	W/S	GU	(91)
	basepair_index[93] = 29	#	trans	W/S	UA	(92)
	basepair_index[94] = 29	#	trans	W/S	UC	(93)
	basepair_index[95] = 30	#	trans	W/S	UG	(94)
	basepair_index[96] = 29	#	trans	W/S	UU	(95)

	#	<*****matrix_7*****>
	basepair_index[97] = 0		#	cis	H/W	AA	(96)
	basepair_index[98] = 0		#	cis	H/W	AC	(97)
	basepair_index[99] = 33	#	cis	H/W	AG	(98)
	basepair_index[100] = 31	#	cis	H/W	AU	(99)
	basepair_index[101] = 0	#	cis	H/W	CA	(100)
	basepair_index[102] = 32	#	cis	H/W	CC	(101)
	basepair_index[103] = 0	#	cis	H/W	CG	(102)
	basepair_index[104] = 0	#	cis	H/W	CU	(103)
	basepair_index[105] = 33	#	cis	H/W	GA	(104)
	basepair_index[106] = 31	#	cis	H/W	GC	(105)
	basepair_index[107] = 34	#	cis	H/W	GG	(106)
	basepair_index[108] = 31	#	cis	H/W	GU	(107)
	basepair_index[109] = 33	#	cis	H/W	UA	(108)
	basepair_index[110] = 31	#	cis	H/W	UC	(109)
	basepair_index[111] = 0	#	cis	H/W	UG	(110)
	basepair_index[112] = 32	#	cis	H/W	UU	(111)

	#	<*****matrix_8*****>
	basepair_index[113] = 38	#	trans	H/W	AA	(112)
	basepair_index[114] = 36	#	trans	H/W	AC	(113)
	basepair_index[115] = 0	#	trans	H/W	AG	(114)
	basepair_index[116] = 35	#	trans	H/W	AU	(115)
	basepair_index[117] = 0	#	trans	H/W	CA	(116)
	basepair_index[118] = 35	#	trans	H/W	CC	(117)
	basepair_index[119] = 0	#	trans	H/W	CG	(118)
	basepair_index[120] = 0	#	trans	H/W	CU	(119)
	basepair_index[121] = 38	#	trans	H/W	GA	(120)
	basepair_index[122] = 36	#	trans	H/W	GC	(121)
	basepair_index[123] = 39	#	trans	H/W	GG	(122)
	basepair_index[124] = 37	#	trans	H/W	GU	(123)
	basepair_index[125] = 0	#	trans	H/W	UA	(124)
	basepair_index[126] = 0	#	trans	H/W	UC	(125)
	basepair_index[127] = 38	#	trans	H/W	UG	(126)
	basepair_index[128] = 36	#	trans	H/W	UU	(127)

	#	<*****matrix_9*****>
	basepair_index[129] = 0	#	cis	H/H	AA	(128)
	basepair_index[130] = 0	#	cis	H/H	AC	(129)
	basepair_index[131] = 41	#	cis	H/H	AG	(130)
	basepair_index[132] = 0	#	cis	H/H	AU	(131)
	basepair_index[133] = 0	#	cis	H/H	CA	(132)
	basepair_index[134] = 0	#	cis	H/H	CC	(133)
	basepair_index[135] = 40	#	cis	H/H	CG	(134)
	basepair_index[136] = 0	#	cis	H/H	CU	(135)
	basepair_index[137] = 41	#	cis	H/H	GA	(136)
	basepair_index[138] = 40	#	cis	H/H	GC	(137)
	basepair_index[139] = 40	#	cis	H/H	GG	(138)
	basepair_index[140] = 0	#	cis	H/H	GU	(139)
	basepair_index[141] = 0	#	cis	H/H	UA	(140)
	basepair_index[142] = 0	#	cis	H/H	UC	(141)
	basepair_index[143] = 0	#	cis	H/H	UG	(142)
	basepair_index[144] = 0	#	cis	H/H	UU	(143)

	#	<*****matrix_10*****>
	basepair_index[145] = 42	#	trans	H/H	AA	(144)
	basepair_index[146] = 42	#	trans	H/H	AC	(145)
	basepair_index[147] = 43	#	trans	H/H	AG	(146)
	basepair_index[148] = 43	#	trans	H/H	AU	(147)
	basepair_index[149] = 42	#	trans	H/H	CA	(148)
	basepair_index[150] = 0	#	trans	H/H	CC	(149)
	basepair_index[151] = 42	#	trans	H/H	CG	(150)
	basepair_index[152] = 43	#	trans	H/H	CU	(151)
	basepair_index[153] = 43	#	trans	H/H	GA	(152)
	basepair_index[154] = 42	#	trans	H/H	GC	(153)
	basepair_index[155] = 44	#	trans	H/H	GG	(154)
	basepair_index[156] = 0	#	trans	H/H	GU	(155)
	basepair_index[157] = 43	#	trans	H/H	UA	(156)
	basepair_index[158] = 43	#	trans	H/H	UC	(157)
	basepair_index[159] = 0	#	trans	H/H	UG	(158)
	basepair_index[160] = 0	#	trans	H/H	UU	(159)

	#	<*****matrix_11*****>
	basepair_index[161] = 45	#	cis	H/S	AA	(160)
	basepair_index[162] = 45	#	cis	H/S	AC	(161)
	basepair_index[163] = 45	#	cis	H/S	AG	(162)
	basepair_index[164] = 45	#	cis	H/S	AU	(163)
	basepair_index[165] = 45	#	cis	H/S	CA	(164)
	basepair_index[166] = 45	#	cis	H/S	CC	(165)
	basepair_index[167] = 45	#	cis	H/S	CG	(166)
	basepair_index[168] = 45	#	cis	H/S	CU	(167)
	basepair_index[169] = 45	#	cis	H/S	GA	(168)
	basepair_index[170] = 0	#	cis	H/S	GC	(169)
	basepair_index[171] = 45	#	cis	H/S	GG	(170)
	basepair_index[172] = 0	#	cis	H/S	GU	(171)
	basepair_index[173] = 46	#	cis	H/S	UA	(172)
	basepair_index[174] = 45	#	cis	H/S	UC	(173)
	basepair_index[175] = 45	#	cis	H/S	UG	(174)
	basepair_index[176] = 45	#	cis	H/S	UU	(175)

	#	<*****matrix_12*****>
	basepair_index[177] = 47	#	trans	H/S	AA	(176)
	basepair_index[178] = 47	#	trans	H/S	AC	(177)
	basepair_index[179] = 47	#	trans	H/S	AG	(178)
	basepair_index[180] = 47	#	trans	H/S	AU	(179)
	basepair_index[181] = 47	#	trans	H/S	CA	(180)
	basepair_index[182] = 47	#	trans	H/S	CC	(181)
	basepair_index[183] = 0	#	trans	H/S	CG	(182)
	basepair_index[184] = 47	#	trans	H/S	CU	(183)
	basepair_index[185] = 0	#	trans	H/S	GA	(184)
	basepair_index[186] = 0	#	trans	H/S	GC	(185)
	basepair_index[187] = 48	#	trans	H/S	GG	(186)
	basepair_index[188] = 0	#	trans	H/S	GU	(187)
	basepair_index[189] = 48	#	trans	H/S	UA	(188)
	basepair_index[190] = 0	#	trans	H/S	UC	(189)
	basepair_index[191] = 48	#	trans	H/S	UG	(190)
	basepair_index[192] = 0	#	trans	H/S	UU	(191)

	#	<*****matrix_13*****>
	basepair_index[193] = 49	#	cis	S/W	AA	(192)
	basepair_index[194] = 50	#	cis	S/W	AC	(193)
	basepair_index[195] = 51	#	cis	S/W	AG	(194)
	basepair_index[196] = 52	#	cis	S/W	AU	(195)
	basepair_index[197] = 49	#	cis	S/W	CA	(196)
	basepair_index[198] = 50	#	cis	S/W	CC	(197)
	basepair_index[199] = 51	#	cis	S/W	CG	(198)
	basepair_index[200] = 52	#	cis	S/W	CU	(199)
	basepair_index[201] = 49	#	cis	S/W	GA	(200)
	basepair_index[202] = 50	#	cis	S/W	GC	(201)
	basepair_index[203] = 53	#	cis	S/W	GG	(202)
	basepair_index[204] = 52	#	cis	S/W	GU	(203)
	basepair_index[205] = 49	#	cis	S/W	UA	(204)
	basepair_index[206] = 50	#	cis	S/W	UC	(205)
	basepair_index[207] = 51	#	cis	S/W	UG	(206)
	basepair_index[208] = 52	#	cis	S/W	UU	(207)

	#	<*****matrix_14*****>
	basepair_index[209] = 54	#	trans	S/W	AA	(208)
	basepair_index[210] = 54	#	trans	S/W	AC	(209)
	basepair_index[211] = 0	#	trans	S/W	AG	(210)
	basepair_index[212] = 56	#	trans	S/W	AU	(211)
	basepair_index[213] = 54	#	trans	S/W	CA	(212)
	basepair_index[214] = 54	#	trans	S/W	CC	(213)
	basepair_index[215] = 55	#	trans	S/W	CG	(214)
	basepair_index[216] = 56	#	trans	S/W	CU	(215)
	basepair_index[217] = 54	#	trans	S/W	GA	(216)
	basepair_index[218] = 54	#	trans	S/W	GC	(217)
	basepair_index[219] = 0	#	trans	S/W	GG	(218)
	basepair_index[220] = 57	#	trans	S/W	GU	(219)
	basepair_index[221] = 54	#	trans	S/W	UA	(220)
	basepair_index[222] = 54	#	trans	S/W	UC	(221)
	basepair_index[223] = 55	#	trans	S/W	UG	(222)
	basepair_index[224] = 56	#	trans	S/W	UU	(223)

	#	<*****matrix_15*****>
	basepair_index[225] = 58	#	cis	S/H	AA	(224)
	basepair_index[226] = 58	#	cis	S/H	AC	(225)
	basepair_index[227] = 58	#	cis	S/H	AG	(226)
	basepair_index[228] = 59	#	cis	S/H	AU	(227)
	basepair_index[229] = 58	#	cis	S/H	CA	(228)
	basepair_index[230] = 58	#	cis	S/H	CC	(229)
	basepair_index[231] = 0	#	cis	S/H	CG	(230)
	basepair_index[232] = 58	#	cis	S/H	CU	(231)
	basepair_index[233] = 58	#	cis	S/H	GA	(232)
	basepair_index[234] = 58	#	cis	S/H	GC	(233)
	basepair_index[235] = 58	#	cis	S/H	GG	(234)
	basepair_index[236] = 58	#	cis	S/H	GU	(235)
	basepair_index[237] = 58	#	cis	S/H	UA	(236)
	basepair_index[238] = 58	#	cis	S/H	UC	(237)
	basepair_index[239] = 0	#	cis	S/H	UG	(238)
	basepair_index[240] = 58	#	cis	S/H	UU	(239)

	#	<*****matrix_16*****>
	basepair_index[241] = 60	#	trans	S/H	AA	(240)
	basepair_index[242] = 60	#	trans	S/H	AC	(241)
	basepair_index[243] = 0	#	trans	S/H	AG	(242)
	basepair_index[244] = 61	#	trans	S/H	AU	(243)
	basepair_index[245] = 60	#	trans	S/H	CA	(244)
	basepair_index[246] = 60	#	trans	S/H	CC	(245)
	basepair_index[247] = 0	#	trans	S/H	CG	(246)
	basepair_index[248] = 0	#	trans	S/H	CU	(247)
	basepair_index[249] = 60	#	trans	S/H	GA	(248)
	basepair_index[250] = 0	#	trans	S/H	GC	(249)
	basepair_index[251] = 61	#	trans	S/H	GG	(250)
	basepair_index[252] = 61	#	trans	S/H	GU	(251)
	basepair_index[253] = 60	#	trans	S/H	UA	(252)
	basepair_index[254] = 60	#	trans	S/H	UC	(253)
	basepair_index[255] = 0	#	trans	S/H	UG	(254)
	basepair_index[256] = 0	#	trans	S/H	UU	(255)

	#	<*****matrix_17*****>
	basepair_index[257] = 62	#	cis	S/S	AA	(256)
	basepair_index[258] = 62	#	cis	S/S	AC	(257)
	basepair_index[259] = 62	#	cis	S/S	AG	(258)
	basepair_index[260] = 62	#	cis	S/S	AU	(259)
	basepair_index[261] = 62	#	cis	S/S	CA	(260)
	basepair_index[262] = 62	#	cis	S/S	CC	(261)
	basepair_index[263] = 62	#	cis	S/S	CG	(262)
	basepair_index[264] = 62	#	cis	S/S	CU	(263)
	basepair_index[265] = 62	#	cis	S/S	GA	(264)
	basepair_index[266] = 62	#	cis	S/S	GC	(265)
	basepair_index[267] = 62	#	cis	S/S	GG	(266)
	basepair_index[268] = 62	#	cis	S/S	GU	(267)
	basepair_index[269] = 62	#	cis	S/S	UA	(268)
	basepair_index[270] = 62	#	cis	S/S	UC	(269)
	basepair_index[271] = 62	#	cis	S/S	UG	(270)
	basepair_index[272] = 62	#	cis	S/S	UU	(271)

	#	<*****matrix_18*****>
	basepair_index[273] = 63	#	trans	S/S	AA	(272)
	basepair_index[274] = 63	#	trans	S/S	AC	(273)
	basepair_index[275] = 63	#	trans	S/S	AG	(274)
	basepair_index[276] = 63	#	trans	S/S	AU	(275)
	basepair_index[277] = 0	#	trans	S/S	CA	(276)
	basepair_index[278] = 0	#	trans	S/S	CC	(277)
	basepair_index[279] = 0	#	trans	S/S	CG	(278)
	basepair_index[280] = 0	#	trans	S/S	CU	(279)
	basepair_index[281] = 64	#	trans	S/S	GA	(280)
	basepair_index[282] = 64	#	trans	S/S	GC	(281)
	basepair_index[283] = 64	#	trans	S/S	GG	(282)
	basepair_index[284] = 64	#	trans	S/S	GU	(283)
	basepair_index[285] = 0	#	trans	S/S	UA	(284)
	basepair_index[286] = 0	#	trans	S/S	UC	(285)
	basepair_index[287] = 0	#	trans	S/S	UG	(286)
	basepair_index[288] = 0	#	trans	S/S	UU	(287)

	return basepair_index


def hash_edge(edge):
	edge = edge.upper()
	if edge == 'W':
		return 0
	elif edge == 'H':
		return 1
	elif edge == 'S':
		return 2
	else:
		return 99999

def reverse_hash_edge(value):
	if value == 0:
		return 'W'
	elif value == 1:
		return 'H'
	elif value == 2:
		return 'S'
	else:
		return ' '

def hash_nucleotide(nucl):
	nucl = nucl.upper()
	if nucl == 'A':
		return 0
	elif nucl == 'C':
		return 1
	elif nucl == 'G':
		return 2
	elif nucl == 'U' or nucl == 'T':
		return 3
	elif nucl == '-':
		return 4
	else:
		return 99999

def reverse_hash_nucleotide(value):
	if value == 0:
		return 'A'
	elif value == 1:
		return 'C'
	elif value == 2:
		return 'G'
	elif value == 3:
		return 'U'
	elif value == 4:
		return '-'
	else:
		return 'N'

def hash_basepair(interaction, ntd_pair):
	edge_index_1 = hash_edge(interaction[1])
	edge_index_2 = hash_edge(interaction[2])
	nuc_index_1 = hash_nucleotide(ntd_pair[0])
	nuc_index_2 = hash_nucleotide(ntd_pair[1])

	is_cis = True
	if interaction[0].lower() == 'c':
		is_cis = True
	elif interaction[0].lower() == 't':
		is_cis = False


	hash_index = (edge_index_1 * 3 + edge_index_2) * 2;
	if is_cis == False:
		hash_index += 1
	hash_index = hash_index * 16 + (nuc_index_1 * 4 + nuc_index_2)

	return hash_index + 1;

def hash_stacking(interaction, ntd_pair):
	if interaction.lower() == 'upward':
		return -1
	elif interaction.lower() == 'downward':
		return -2
	elif interaction.lower() == 'inward':
		return -3
	elif interaction.lower() == 'outward':
		return -4
	else:
		return 99999

def reverse_hash_stacking(value):
	if value == -1:
		return 'upward'
	elif value == -2:
		return 'downward'
	elif value == -3:
		return 'inward'
	elif value == -4:
		return 'outward'
	else:
		return ' '

def interpret_hash_value(value):

	# value = int(value)
	if value < 0:
		# base-stacking
		if value == -1:
			return 'upward', ''
		elif value == -2:
			return 'downward', ''
		elif value == -3:
			return 'inward', ''
		elif value == -4:
			return 'outward', ''
		
	elif value == 0:
		# hydrogen bond
		return 'hbond', ''
	else:
		# base-pair
	
		edge_and_orientation = math.floor((value - 1) / 16)

		ntd_pair_value = value - 1 - edge_and_orientation * 16
		nuc_index_1 = math.floor(ntd_pair_value / 4)
		nuc_index_2 = ntd_pair_value % 4

		ntd_pair = reverse_hash_nucleotide(nuc_index_1) + reverse_hash_nucleotide(nuc_index_2)
		
		is_cis = True
		if edge_and_orientation % 2 == 1:
			is_cis = False

		edges = math.floor(edge_and_orientation / 2)
		edge_1 = math.floor(edges / 3)
		edge_2 = edges % 3

		interaction = reverse_hash_edge(edge_1) + reverse_hash_edge(edge_2)

		if is_cis:
			interaction = 'c' + interaction
		else:
			interaction = 't' + interaction

		return interaction, ntd_pair

	return 'N/A', 'N/A'
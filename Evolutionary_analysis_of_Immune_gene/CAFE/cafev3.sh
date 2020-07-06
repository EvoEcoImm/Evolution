~/opt/CAFE/release/cafe
load -i allmanrefine_count_cafe.tab -t 8 -l modelchoose -p 0.05
tree -i cafetree.tre
esterror -dataerror allmanrefine_errorzntrans_count_cafe.csv allmanrefine_errorzng_count_cafe.csv -diff 9 -o Error_znboth_all.txt
errormodel -model Error_znboth_all.txt -all

#k-means test
(1,(1,((1,1)1,(1,(1,((1,(1,1)1)1,(1,((1,1)1,(1,(1,(1,((1,1)1,1)1)1)1)1)1)1)1)1)1)1)1)
lambda -s -t (2,(1,((1,1)1,(1,(1,((1,(1,1)1)1,(1,((1,1)1,(1,(1,(1,((1,1)1,1)1)1)1)1)1)1)1)1)1)1)1)
lambda -s -t (1,(2,((1,1)1,(1,(1,((1,(1,1)1)1,(1,((1,1)1,(1,(1,(1,((1,1)1,1)1)1)1)1)1)1)1)1)1)1)1)
lambda -s -t (1,(1,((2,1)1,(1,(1,((1,(1,1)1)1,(1,((1,1)1,(1,(1,(1,((1,1)1,1)1)1)1)1)1)1)1)1)1)1)1)
lambda -s -t (1,(1,((1,2)1,(1,(1,((1,(1,1)1)1,(1,((1,1)1,(1,(1,(1,((1,1)1,1)1)1)1)1)1)1)1)1)1)1)1)
lambda -s -t (1,(1,((1,1)2,(1,(1,((1,(1,1)1)1,(1,((1,1)1,(1,(1,(1,((1,1)1,1)1)1)1)1)1)1)1)1)1)1)1)
lambda -s -t (1,(1,((1,1)1,(2,(1,((1,(1,1)1)1,(1,((1,1)1,(1,(1,(1,((1,1)1,1)1)1)1)1)1)1)1)1)1)1)1)
lambda -s -t (1,(1,((1,1)1,(1,(2,((1,(1,1)1)1,(1,((1,1)1,(1,(1,(1,((1,1)1,1)1)1)1)1)1)1)1)1)1)1)1)
lambda -s -t (1,(1,((1,1)1,(1,(1,((2,(1,1)1)1,(1,((1,1)1,(1,(1,(1,((1,1)1,1)1)1)1)1)1)1)1)1)1)1)1)
lambda -s -t (1,(1,((1,1)1,(1,(1,((1,(2,1)1)1,(1,((1,1)1,(1,(1,(1,((1,1)1,1)1)1)1)1)1)1)1)1)1)1)1)
lambda -s -t (1,(1,((1,1)1,(1,(1,((1,(1,2)1)1,(1,((1,1)1,(1,(1,(1,((1,1)1,1)1)1)1)1)1)1)1)1)1)1)1)
lambda -s -t (1,(1,((1,1)1,(1,(1,((1,(1,1)2)1,(1,((1,1)1,(1,(1,(1,((1,1)1,1)1)1)1)1)1)1)1)1)1)1)1)
lambda -s -t (1,(1,((1,1)1,(1,(1,((1,(1,1)1)2,(1,((1,1)1,(1,(1,(1,((1,1)1,1)1)1)1)1)1)1)1)1)1)1)1)
lambda -s -t (1,(1,((1,1)1,(1,(1,((1,(1,1)1)1,(2,((1,1)1,(1,(1,(1,((1,1)1,1)1)1)1)1)1)1)1)1)1)1)1)
lambda -s -t (1,(1,((1,1)1,(1,(1,((1,(1,1)1)1,(1,((2,1)1,(1,(1,(1,((1,1)1,1)1)1)1)1)1)1)1)1)1)1)1)
lambda -s -t (1,(1,((1,1)1,(1,(1,((1,(1,1)1)1,(1,((1,2)1,(1,(1,(1,((1,1)1,1)1)1)1)1)1)1)1)1)1)1)1)
lambda -s -t (1,(1,((1,1)1,(1,(1,((1,(1,1)1)1,(1,((1,1)2,(1,(1,(1,((1,1)1,1)1)1)1)1)1)1)1)1)1)1)1)
lambda -s -t (1,(1,((1,1)1,(1,(1,((1,(1,1)1)1,(1,((1,1)1,(2,(1,(1,((1,1)1,1)1)1)1)1)1)1)1)1)1)1)1)
lambda -s -t (1,(1,((1,1)1,(1,(1,((1,(1,1)1)1,(1,((1,1)1,(1,(2,(1,((1,1)1,1)1)1)1)1)1)1)1)1)1)1)1)
lambda -s -t (1,(1,((1,1)1,(1,(1,((1,(1,1)1)1,(1,((1,1)1,(1,(1,(2,((1,1)1,1)1)1)1)1)1)1)1)1)1)1)1)
lambda -s -t (1,(1,((1,1)1,(1,(1,((1,(1,1)1)1,(1,((1,1)1,(1,(1,(1,((2,1)1,1)1)1)1)1)1)1)1)1)1)1)1)
lambda -s -t (1,(1,((1,1)1,(1,(1,((1,(1,1)1)1,(1,((1,1)1,(1,(1,(1,((1,2)1,1)1)1)1)1)1)1)1)1)1)1)1)
lambda -s -t (1,(1,((1,1)1,(1,(1,((1,(1,1)1)1,(1,((1,1)1,(1,(1,(1,((1,1)2,1)1)1)1)1)1)1)1)1)1)1)1)
lambda -s -t (1,(1,((1,1)1,(1,(1,((1,(1,1)1)1,(1,((1,1)1,(1,(1,(1,((1,1)1,2)1)1)1)1)1)1)1)1)1)1)1)
lambda -s -t (1,(1,((1,1)1,(1,(1,((1,(1,1)1)1,(1,((1,1)1,(1,(1,(1,((1,1)1,1)2)1)1)1)1)1)1)1)1)1)1)
lambda -s -t (1,(1,((1,1)1,(1,(1,((1,(1,1)1)1,(1,((1,1)1,(1,(1,(1,((1,1)1,1)1)2)1)1)1)1)1)1)1)1)1)
lambda -s -t (1,(1,((1,1)1,(1,(1,((1,(1,1)1)1,(1,((1,1)1,(1,(1,(1,((1,1)1,1)1)1)2)1)1)1)1)1)1)1)1)
lambda -s -t (1,(1,((1,1)1,(1,(1,((1,(1,1)1)1,(1,((1,1)1,(1,(1,(1,((1,1)1,1)1)1)1)2)1)1)1)1)1)1)1)
lambda -s -t (1,(1,((1,1)1,(1,(1,((1,(1,1)1)1,(1,((1,1)1,(1,(1,(1,((1,1)1,1)1)1)1)1)2)1)1)1)1)1)1)
lambda -s -t (1,(1,((1,1)1,(1,(1,((1,(1,1)1)1,(1,((1,1)1,(1,(1,(1,((1,1)1,1)1)1)1)1)1)2)1)1)1)1)1)
lambda -s -t (1,(1,((1,1)1,(1,(1,((1,(1,1)1)1,(1,((1,1)1,(1,(1,(1,((1,1)1,1)1)1)1)1)1)1)2)1)1)1)1)
lambda -s -t (1,(1,((1,1)1,(1,(1,((1,(1,1)1)1,(1,((1,1)1,(1,(1,(1,((1,1)1,1)1)1)1)1)1)1)1)2)1)1)1)
lambda -s -t (1,(1,((1,1)1,(1,(1,((1,(1,1)1)1,(1,((1,1)1,(1,(1,(1,((1,1)1,1)1)1)1)1)1)1)1)1)2)1)1)
lambda -s -t (1,(1,((1,1)1,(1,(1,((1,(1,1)1)1,(1,((1,1)1,(1,(1,(1,((1,1)1,1)1)1)1)1)1)1)1)1)1)2)1)
lambda -s -t (1,(1,((1,1)1,(1,(1,((1,(1,1)1)1,(1,((1,1)1,(1,(1,(1,((1,1)1,1)1)1)1)1)1)1)1)1)1)1)2)

#model choose
lambda -s -t (1,(1,((2,2)2,(3,(2,((2,(2,2)2)2,(2,((3,3)3,(4,(4,(4,((4,4)4,4)4)4)4)4)3)2)2)2)2)2)1)
lambda -s -t (1,(1,((2,2)2,(3,(2,((2,(2,2)2)2,(2,((3,3)3,(3,(3,(3,((3,3)3,3)3)3)3)3)3)2)2)2)2)2)1)
lambda -s -t (1,(1,((2,2)2,(3,(3,((3,(3,3)3)3,(4,((4,4)4,(5,(5,(5,((5,5)5,5)5)5)5)5)4)4)3)3)3)2)1)
lambda -s -t (1,(1,((2,2)2,(3,(3,((3,(3,3)3)3,(3,((3,3)3,(4,(4,(4,((4,4)4,4)4)4)4)4)3)3)3)3)3)2)1)
lambda -s -t (1,(1,((2,2)2,(3,(3,((3,(3,3)3)3,(3,((3,3)3,(3,(3,(3,((3,3)3,3)3)3)3)3)3)3)3)3)3)2)1)
lambda -s -t (1,(1,((2,2)2,(2,(2,((2,(2,2)2)2,(2,((2,2)2,(2,(2,(2,((2,2)2,2)2)2)2)2)2)2)2)2)2)2)1)
lambda -s -t (1,(1,((2,2)2,(2,(2,((2,(2,2)2)2,(3,((3,3)3,(3,(3,(3,((3,3)3,3)3)3)3)3)2)2)2)2)2)2)1)
lambda -s -t (1,(1,((1,1)1,(2,(2,((2,(2,2)2)2,(2,((2,2)2,(2,(2,(2,((2,2)2,2)2)2)2)2)2)2)2)2)2)1)1)
lambda -s
lambda -s -t (1,(2,((3,4)5,(6,(7,((8,(9,10)11)12,(13,((14,15)16,(17,(18,(19,((20,21)22,23)24)25)26)27)28)29)30)31)32)33)34)
(1,(2,((3,4)5,(6,(7,((8,(9,10)11)12,(13,((14,15)16,(17,(18,(19,((20,21)22,23)24)25)26)27)28)29)30)31)32)33)34)
#k-means
lambda -s -t (1,(1,((1,1)1,(1,(1,((1,(1,1)1)1,(1,((1,1)1,(1,(1,(1,((2,1)2,1)1)1)1)1)1)1)1)2)2)2)1)
lambda -s -t (3,(3,((3,1)3,(3,(3,((3,(3,3)3)3,(3,((3,3)3,(3,(3,(3,((1,3)2,3)3)3)1)1)3)1)3)1)1)1)3)
lambda -s -t (1,(1,((1,2)3,(3,(3,((3,(1,3)3)3,(3,((3,1)3,(3,(3,(3,((2,1)4,3)1)3)1)2)1)1)3)2)2)2)1)
lambda -s -t (1,(1,((1,1)5,(3,(5,((5,(3,5)5)3,(5,((3,1)5,(5,(5,(5,((2,3)4,3)3)5)1)2)1)1)5)2)2)2)3)
lambda -s -t (6,(6,((6,1)5,(3,(5,((5,(3,5)5)3,(5,((3,6)5,(5,(5,(5,((1,3)4,3)3)5)6)2)6)6)5)1)2)1)3)
lambda -s -t (7,(7,((7,1)5,(3,(5,((5,(6,5)5)3,(5,((3,6)5,(5,(5,(5,((1,6)4,3)3)5)7)2)7)7)5)1)2)1)6)
lambda -s -t (7,(7,((7,1)5,(8,(3,((3,(6,5)3)8,(3,((8,6)3,(5,(5,(3,((1,6)4,8)8)3)7)2)7)7)3)1)2)1)6)
lambda -s -t (7,(7,((7,1)5,(8,(3,((3,(6,5)3)8,(3,((8,7)3,(5,(5,(3,((2,6)4,8)6)3)1)9)7)1)3)2)9)2)6)
lambda -s -t (7,(7,((7,1)5,(10,(3,((3,(6,5)3)10,(3,((10,7)3,(8,(8,(3,((2,6)4,10)6)3)1)9)7)1)3)2)9)2)6)
lambda -s -t (7,(7,((7,1)5,(10,(11,((3,(6,5)3)10,(11,((10,7)3,(8,(8,(3,((2,6)4,10)6)3)1)9)7)1)3)2)9)2)6)
lambda -s -t (7,(7,((7,1)11,(10,(5,((3,(6,11)3)10,(5,((10,7)3,(12,(12,(8,((2,6)4,10)6)3)1)9)7)1)3)2)9)2)6)
lambda -s -t (1,(1,((1,2)11,(10,(5,((3,(6,11)3)10,(5,((10,7)3,(12,(12,(8,((9,6)4,10)6)3)2)13)1)2)3)9)13)9)6)
lambda -s -t (1,(1,((1,14)11,(10,(5,((3,(6,11)3)10,(5,((10,7)3,(12,(12,(8,((9,6)4,10)6)3)2)13)1)2)3)9)13)9)6)
lambda -s -t (1,(1,((1,14)15,(10,(5,((3,(6,15)3)10,(5,((10,7)3,(12,(12,(11,((9,6)4,10)6)3)2)13)1)2)3)9)13)9)6)
lambda -s -t (1,(1,((1,14)15,(16,(5,((3,(6,15)3)16,(5,((16,7)3,(12,(12,(11,((9,6)4,16)10)3)2)13)1)2)3)9)13)9)6)
lambda -s -t (1,(1,((1,14)15,(16,(5,((3,(6,15)3)10,(5,((10,7)3,(12,(12,(11,((9,6)4,10)17)3)2)13)1)2)3)9)13)9)6)


#lhest-test
load -i allmanrefine_count_cafe.tab -t 4 -l simcafelog -p 0.05
tree -i cafetree.tre
lambda -l 0.00207147786220
genfamily smgenefamily/rnd -t 100
lhtest -d smgenefamily -l 0.00207147786220 -t (2,(2,((2,2)3,(3,(3,((3,(2,3)3)3,(3,((3,2)3,(3,(3,(3,((1,2)3,3)3)3)2)1)3)2)3)3)1)1)2) -o lhtest_result

#output
load -i allmanrefine_count_cafe.tab -t 4 -l cafelog -p 0.05
tree -i cafetree.tre
errormodel -model Error_znboth_all.txt -all
lambda -s -t (1,(1,((2,2)2,(2,(2,((2,(2,2)2)2,(2,((2,2)2,(2,(2,(2,((2,2)2,2)2)2)2)2)2)2)2)2)2)2)1)
report output_result2

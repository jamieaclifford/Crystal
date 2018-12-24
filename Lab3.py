
from numpy import *
from scipy.signal import correlate2d
from numpy.random import randint,choice,uniform
from matplotlib.pyplot import *
import matplotlib.pyplot as plt

from os import system

class KineticMonteCarlo(object) :

    def __init__( self, Model ) :

        # reference state
        self.Model = Model
        self.state = Model.state

        # event rates and evolution rule
        self.evolve = Model.evolve
        self.get_rates = Model.get_rates

    def time_step( self ) :

        # calculate transition rate fractions
        rates = self.get_rates()  ### a k-dimensional array
        probabilities = cumsum(rates)  #### a k-dim array
	print('Probabilities are'), probabilities
        total_rate = probabilities[-1]
        probabilities /= total_rate

        # choose events according to rates
        u = uniform(0,1)
        #get array of true and false depending on value of random number which is compared to the intervals for events
        events_index = sum(probabilities < u)
	if events_index==0 :
		print ('Event is adsorption')
	if events_index==1 :
		print ('Event is diffusion')
	if events_index==2 :
		print ('Event is agglomeration')

        # generate waiting time
        v = uniform(0,1)
        dt = -log(v)/total_rate
        # carry out events
        self.evolve(events_index,dt)

class LatticeGas(object) :

    def __init__(self, n_sites, density ) :

        # initialisation of lattice gas with given denisty
        self.size = 2*int(ceil(sqrt(n_sites)/2))
        self.state = choice([0,1],size=(self.size,self.size),p=[1-density,density])
        self.time = 0.
        self.n=0
        self.singles=0
        self.data=[]
        self.atom=[]

        self.zero=[]

        self.sing=[]

    def alone(self):

        self.n +=1
        array=()
        self.state==1
        mask = self.state == 1
        ilist,jlist = where (mask)
        oneslist = zip (ilist,jlist)
        k=randint(len(oneslist))
        oneslist [k]
        print oneslist [k]
        row=ilist [k]
        column=jlist [k]


        nearn=self.state[(row - 1)%self.size,column%self.size]+self.state[(row + 1)%self.size,column%self.size]+self.state[row%self.size,(column-1)%self.size]+self.state[row%self.size,(column+1)%self.size]
        if nearn==0:
            self.singles=self.singles+1
        else:
            self.singles=self.singles+0

        return self.singles

    def empty(self):
        self.n +=1
        array1=()
        self.state==0
        mask = self.state == 0
        ilist,jlist = where (mask)
        zerolist= zip (ilist,jlist)
        self.number_emp=len(zerolist)
        print('0=',self.number_emp)
        return self.number_emp

    def islands(self) :

        #Helper function that indicates if we can enter the cell or not
        def can_enter_cell(matrix, is_visited, cur_row, cur_col) :
            n_rows = len(matrix)
            n_cols = len(matrix[0])

            #If we are outside the bounds of the matrix or
            #if the cell is already visited or if the value in cell is 0
            #then we shouldn't enter the cell
            if (cur_row < 0 or cur_row >= n_rows
                or cur_col < 0 or cur_col >= n_cols
                or is_visited[cur_row][cur_col]
                or matrix[cur_row][cur_col] == 0) :
                return False
            return True

        #Helper function to count the number of islands of 1's
        #matrix: 2-D matrix consisting of 0's and 1's
        #is_visited: if cell (i, j) has been visited, is_visited[i][j] is set to True
        #cur_row: row of the current cell being processed
        #cur_col: column of the current cell being processed
        def expand_search(matrix, is_visited, cur_row, cur_col) :
            n_rows = len(matrix)
            n_cols = len(matrix[0])

            #at the point a cell is visited via expand it will be referred to as vsited
            is_visited[cur_row][cur_col] = True

            #For the current cell, find out if we can continue the island of 1's
            #with its neighbors. Each cell has 8 neighbors. The rows
            #of neighbors will vary from cur_row - 1 to cur_row + 1
            #The columns of the neighbors will vary from cur_col - 1
            #to cur_col + 1
            for  i in range(-1, 2):
                for  j in range(-1, 2):
                    #cell part of island
                    is_safe_cell = can_enter_cell(matrix, is_visited, cur_row+i,
                                        cur_col+j)

                    if (is_safe_cell) :
                        #check the neighbours for all the ranges up and down
                        expand_search(matrix, is_visited, cur_row+i, cur_col+j)
        #Main function to find the number of islands of 1's
        #matrix: 2-D matrix consisting of 0's and 1's. Should not be empty
        def find_islands(matrix) :
            n_rows = len(matrix)
            n_cols = len(matrix[0])
            #list of all the points visited as true or false (boolean)
            is_visited = [ [False for x in range(n_cols)] for x in range(n_rows)]

            #Search all the cells in matrix that are not yet visited
            count = 0
            for  i in range(0, n_rows):
                for  j in range(0, n_cols):
                    if (matrix[i][j] == 1 and not is_visited[i][j]) :
                        #We assume we have found an island and use the atom as starting point. Now expand the island
                        #in all directions
                        #look through all the points in the lattice(represented by a matrix)
                        island = 1
                        def alone():
                            for i_alone in range(-1, 1):
                                for j_alone in range(-1, 1):
                                    if (i_alone < 0 or i_alone >= len(matrix)
                                        or j_alone < 0 or j_alone >= len(matrix[0])):
                                            if(matrix[i_alone][j_alone]==1):
                                                return False
                            return True
                        if(not alone()):
                            island = 1
                            expand_search(matrix, is_visited, i, j)
                        else:
                            island = 0
                        count = count + island
            return count

        islands = 0
        tmp = self.state
        array=()
        tmp==1
        mask = tmp == 1
        ilist,jlist = where (mask)
        oneslist = zip (ilist,jlist)
        tmpMap = []
        for i in range(32): #change for sqrt of size
            tmpRow = []
            for j in range(32):
                #fill temp matrix with zeros
                tmpRow.append(0)
            tmpMap.append(tmpRow)
        i = 0
        for one in oneslist:
            #go through the duplicated lattice and access all coordinates
            tmpMap[one[0]][one[1]] = 1
        #print(tmpMap)
        return find_islands(tmpMap)

    def get_rates(self) :

        prefactor_ad=0.001
        energybarrier_ad=0.01

        prefactor_d=0.1
        energybarrier_d=0.35

        prefactor_ag=0.1
        energybarrier_ag=0.5


        self.rate_ad=prefactor_ad*exp(-energybarrier_ad)
        self.rate_d=prefactor_d*exp(-energybarrier_d)
        self.rate_ag=prefactor_ag*exp(-energybarrier_ag)

        #assuming B=1.0

        return array([self.rate_ad,self.rate_d,self.rate_ag])

    def evolve(self,events_index,dt) :
		self.n += 1
		if events_index == 0 :
			#adsorption
			n=self.size

			row_a = randint(0,self.size)
			column_a = randint(0,self.size)

			print row_a
			print column_a
			print self.state[row_a,column_a]
			if self.state[row_a,column_a] == 0:
				self.state[row_a,column_a]=1
			print self.state

			self.time= self.time + dt/self.size**2
			print self.time

		if events_index == 1 :
			#diffusion
			array=()
			self.state==1
			mask = self.state == 1
			ilist,jlist = where (mask)
			oneslist = zip (ilist,jlist)
			k=randint(len(oneslist))
			oneslist [k]
			print oneslist [k]

			row=ilist [k]
			column=jlist [k]

			nearn=self.state[(row - 1)%self.size,column%self.size]+self.state[(row + 1)%self.size,column%self.size]+self.state[row%self.size,(column-1)%self.size]+self.state[row%self.size,(column+1)%self.size]

			if nearn > 0 :
				print ("Diffusion not possible")

			if nearn ==0 :
				print ("Diffusion is possible")
				move= []
				#Checking nearest neighbours of up
				nearn_up=self.state[(row- 2)%self.size,column%self.size]+self.state[(row-1)%self.size,(column-1)%self.size]+self.state[(row-1)%self.size,(column+1)%self.size]
				if nearn_up > 0:
					print("Cannot diffuse up")
				if nearn_up == 0 :
					move.append('up')
					print ("Diffusion possible up")

				#Checking nearet neighbours of right
				nearn_right=self.state[(row- 1)%self.size,(column+1)%self.size]+self.state[row%self.size,(column+2)%self.size]+self.state[(row+1)%self.size,(column+1)%self.size]
				if nearn_right > 0:
					print("Cannot diffuse right")
				if nearn_right == 0 :
					move.append('right')
					print ("Diffusion possible right")

				#Checking nearest neighbours of down
				nearn_down=self.state[(row+ 2)%self.size,(column-1)%self.size]+self.state[(row+2)%self.size,column%self.size]+self.state[(row+1)%self.size,(column+1)%self.size]
				if nearn_down > 0:
					print("Cannot diffuse down")
				if nearn_down == 0 :
					move.append('down')
					print ("Diffusion possible down")

				#Checking nearest neighbours of left
				nearn_left=self.state[(row- 1)%self.size,(column-1)%self.size]+self.state[row%self.size,(column-2)%self.size]+self.state[(row+1)%self.size,(column-1)%self.size]
				if nearn_left > 0:
					print("Cannot diffuse left")
				if nearn_left == 0 :
					move.append('left')
					print ("Diffusion possible left")

				print move
				if(len(move)>0):
					chosen_move=choice(move)
					print ('Chosen move is:',chosen_move)
					if chosen_move=='up' :
						self.state[row%self.size,column%self.size]=0
						self.state[(row-1)%self.size,column%self.size]=1

					elif chosen_move=='right' :
						self.state[row%self.size,column%self.size]=0
						self.state[row%self.size,(column+1)%self.size]=1

					elif chosen_move=='down' :
						self.state[row%self.size,column%self.size]=0
						self.state[(row +1)%self.size,column%self.size]=1

					elif chosen_move=='left' :
						self.state[row%self.size,column%self.size]=0
						self.state[row%self.size,(column-1)%self.size]=1
				else:
					self.state[row%self.size,column%self.size]=1
					print ('cannot move')

				self.time= self.time + dt/self.size**2
				print self.time

				return self.state

		if events_index == 2 :

			#agglomeration
			array=()
			self.state==1
			mask = self.state == 1
			ilist,jlist = where (mask)
			oneslist = zip (ilist,jlist)
			k=randint(len(oneslist))
			oneslist [k]
			print oneslist [k]

			row=ilist [k]
			column=jlist [k]

			nearn=self.state[(row - 1)%self.size,column%self.size]+self.state[(row + 1)%self.size,column%self.size]+self.state[row%self.size,(column-1)%self.size]+self.state[row%self.size,(column+1)%self.size]
			if nearn > 0 :
				print ("Agglomeration not possible")

			if nearn ==0 :
				print ("Agglomeration is possible")

				move=[]
				#Checking nearest neighbours of up
				nearn_up=self.state[(row- 2)%self.size,column%self.size]+self.state[(row-1)%self.size,(column-1)%self.size]+self.state[(row-1)%self.size,(column+1)%self.size]
				if nearn_up > 0:
					print("Agglomeration posible up")
					move.append('up')
				if nearn_up == 0 :
					print ("Cannot Agglomerate up")

				#Checking nearet neighbours of right
				nearn_right=self.state[(row- 1)%self.size,(column+1)%self.size]+self.state[row%self.size,(column+2)%self.size]+self.state[(row+1)%self.size,(column+1)%self.size]
				if nearn_right > 0:
					print("Agglomeration posible right")
					move.append('right')
				if nearn_right == 0 :
					print ("Cannot Agglomerate right")

				#Checking nearest neighbours of down
				nearn_down=self.state[(row+ 2)%self.size,(column-1)%self.size]+self.state[(row+2)%self.size,column%self.size]+self.state[(row+1)%self.size,(column+1)%self.size]
				if nearn_down > 0:
					print("Agglomeration posible down")
					move.append('down')
				if nearn_down == 0 :
					print ("Cannot Agglomerate down")

				#Checking nearest neighbours of left
				nearn_left=self.state[(row- 1)%self.size,(column-1)%self.size]+self.state[row%self.size,(column-2)%self.size]+self.state[(row+1)%self.size,(column-1)%self.size]
				if nearn_left > 0:
					print("Agglomeration posible left")
					move.append('left')
				if nearn_left == 0 :
					print ("Cannot Agglomerate left")


				print move
				if(len(move)>0):
					chosen_move=choice(move)
					print ('Chosen move is:',chosen_move)

					if chosen_move=='up' :
						self.state[row%self.size,column%self.size]=0
						self.state[(row-1)%self.size,column%self.size]=1

					elif chosen_move=='right' :
						self.state[row%self.size,column%self.size]=0
						self.state[row%self.size,(column+1)%self.size]=1

					elif chosen_move=='down' :
						self.state[row%self.size,column%self.size]=0
						self.state[(row +1)%self.size,column%self.size]=1

					elif chosen_move=='left' :
						self.state[row%self.size,column%self.size]=0
						self.state[row%self.size,(column-1)%self.size]=1

				else:
					self.state[row%self.size,column%self.size]=1
					print("cannot move")

				self.time= self.time + dt/self.size**2
				print self.time
    def atom(self):
        nearn=self.state[(row - 1)%self.size,column%self.size]+self.state[(row + 1)%self.size,column%self.size]+self.state[row%self.size,(column-1)%self.size]+self.state[row%self.size,(column+1)%self.size]
        if nearn==0 :
            self.atom.append((self.time,str(nearn)))
            print self.atom

    def show(self) :

        figure(figsize=(5,5))
        suptitle('Time={:.2e}'.format(self.time))
        title('#Islands(t)= '+str(model.islands()))
        self.data.append((self.time,model.islands()))
        self.sing.append((self.time,model.alone()))
        self.zero.append((self.time,model.empty()))

        imshow(self.state,interpolation="none",cmap='spring_r',vmin=0,vmax=1)
        savefig(str(self.n).zfill(5)+'.png')
        close()


system('rm animation.gif')
# example of how code should work
model = LatticeGas(1024,0.3)
print model.state
#model.show()
kmc = KineticMonteCarlo(model)

for n in range(0,100) :
    	kmc.time_step()
    	model.show()

#plot of empty empty space against time

f=plt.figure()
x, y = zip(*model.data)
plt.scatter(x, y)
plt.savefig('Islands_vs_Time.jpg')
f.show()

g=plt.figure()
x,y=zip(*model.sing)
plt.scatter(x,y)
plt.savefig('single_vs_time.jpg')
g.show()

h=plt.figure()
x,y=zip(*model.zero)
plt.plot(x,y)
plt.savefig('EmptySpace_vs_time.jpg')
h.show()






system('convert -delay 2 -loop 1 *.png animation.gif')
system('rm *.png')

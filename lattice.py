import matplotlib.pyplot as plt
import matplotlib.patches as patches
from math import *
import numpy as np

class Lattice():
    """
    Lattice class for 2D integer lattices.
    Currently supports only 2D lattices.
    """

    #v1 and v2 must be orthogonal
    #WLOG we impose that v1 and v2 have argument in the interval (-pi/2,pi/2] (This is equivalent to imposing the first coordinate be positive, and this will be how we normalize inputs for any dimension)
    #also note that then necessarily if we had n vectors in n dim'l space that were linearly independent all with a coordiante equal to zero, then we get a contradtion with linear dependendance in n-1 dim'l space.
    #Hence necessarily there exists a vector with POSITIVE first coordinate

    def __init__(self,v1,v2,coord=0):
        self.vi = [v1,v2] # to make more general, have the input be vi rather than v1,v2
        self.coord = coord # this is the coordinate we force to have non-negative values for our calculations
        self.make_coordinate_positive() # normalizes the input
        # assume values that we want positive are positive
        self.abs_v1 = sqrt(sum([(v1[i])**2 for i in range(len(v1))]))
        self.abs_v2 = sqrt(sum([(v2[i])**2 for i in range(len(v2))]))
        
    def return_coordinates_of_fundamental_parallelogram(self):
        #ONLY WORKS FOR 2D CASE
        v1 = self.vi[0]
        v2 = self.vi[1]
        average = self.av1_plus_bv2(0.5,0.5,v1,v2)
        x_values = [average[0]+self.av1_plus_bv2(-0.5,0.5,v1,v2)[0],average[0]+self.av1_plus_bv2(-0.5,-0.5,v1,v2)[0],average[0]+self.av1_plus_bv2(0.5,-0.5,v1,v2)[0],average[0]+self.av1_plus_bv2(0.5,0.5,v1,v2)[0]]
        y_values = [average[1]+self.av1_plus_bv2(-0.5,0.5,v1,v2)[1],average[1]+self.av1_plus_bv2(-0.5,-0.5,v1,v2)[1],average[1]+self.av1_plus_bv2(0.5,-0.5,v1,v2)[1],average[1]+self.av1_plus_bv2(0.5,0.5,v1,v2)[1]]
        return x_values, y_values
    
    def return_area_of_fundamental_parallelogram(self):
        v1 = self.vi[0]
        v2 = self.vi[1]
        determinant = v1[0]*v2[1] - v1[1]*v2[0]
        return abs(determinant)


    def find_index_of_vector_with_coord_non_zero(self):
        for i in range(len(self.vi)):
            if self.vi[i][self.coord] > 0:
                return i
            
    def return_vector_with_minimal_non_zero_size_for_coord_i(self,coord_i):
        j = None
        min_size = None
        for i in range(len(self.vi)):
            size = abs(self.vi[i][coord_i])
            if ((min_size == None) or (min_size > size)) and size > 0:
                j = i
                min_size = size
        return self.vi[j], min_size

    def size_of_vector(self,v):
        return sqrt(sum((element**2 for element in v)))
    
    def square_of_size_of_vector(self,v):
        return sum((element**2 for element in v))
    
    def av1_plus_bv2(self,a,b,v1,v2):
        v3 = []
        for i in range(len(v1)):
            v3.append(a*v1[i]+b*v2[i])
        return v3

    def multiply_vector_by_scalar(self,v1,scalar):
        v2 = []
        for i in v1:
            v2.append(scalar*i)
        return v2

    def make_coordinate_positive(self):
        for i in range(len(self.vi)):
            if self.vi[i][self.coord] < 0:
                self.vi[i] = self.multiply_vector_by_scalar(self.vi[i],-1)

    def lattice_points_contained_by_circle_of_radius_R(self,R,bp = True):
        # ONLY WORKS FOR 2D
        interior_points = []
        boundary_points = []
        integer_R = floor(R)
        v1, v1_coord = self.return_vector_with_minimal_non_zero_size_for_coord_i(self.coord)
        for i in range(floor(-(integer_R+1)/v1_coord),ceil((integer_R+1)/v1_coord)):
            v2, v2_coord = self.return_vector_with_minimal_non_zero_size_for_coord_i((self.coord+1)%2)
            num_iterations = 2*ceil(pi*(R**2)/(self.return_area_of_fundamental_parallelogram()))
            for j in range(-num_iterations,num_iterations):
                size = self.size_of_vector(self.av1_plus_bv2(i,j,self.vi[0],self.vi[1])) # this checks size with v_1 and v_2, shoudl check size with fund para vectors
                if size < R:
                    interior_points.append([i,j])
                    print(i,j)
                if size == R:
                    if not bp:
                        interior_points.append([i,j])
                        continue
                    boundary_points.append([i,j])
        if not bp:
            return interior_points
        return interior_points, boundary_points
    
    def compute_points_from_coefficients(self,coefficients):
        #ONLY WORKS FOR 2D ATM
        points = []
        for coefficient_vector in coefficients:
            points.append(self.av1_plus_bv2(coefficient_vector[0],coefficient_vector[1],self.vi[0],self.vi[1]))
        return points
    
    def plot_lattice_points_contained_by_circle_of_radius_R(self,radius,j=0,D=1,tiled=False):
        # tiled = True is just if you're interested in seeing a visual highlighting the periodicity of points with respect to D
        # shouldnt use this function for large R compared to vector sizes
        interior_coefficients, boundary_coefficients = self.lattice_points_contained_by_circle_of_radius_R(radius)
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

        # Major ticks every 10, minor ticks every 1
        major_ticks = np.arange(-2*radius - (-2*radius)%10, 2*radius -(2*radius)%10 +11, 10)
        minor_ticks = np.arange(-2*radius - (-2*radius)%10, 2*radius -(2*radius)%10 +11, 1)

        ax.set_xticks(major_ticks)
        ax.set_xticks(minor_ticks, minor=True)
        ax.set_yticks(major_ticks)
        ax.set_yticks(minor_ticks, minor=True)

        # And a corresponding grid
        ax.grid(which='both')

        # Or if you want different settings for the grids:
        ax.grid(which='minor', alpha=0.2)
        ax.grid(which='major', alpha=0.5)

        angle = [2*pi/100*i for i in range(101)]

        x = [radius*cos(i) for i in angle]
        y = [radius*sin(i) for i in angle]

        plt.axis("equal")
        plt.plot(x, y, color='black')
        print(interior_coefficients)
        i_points,b_points = self.compute_points_from_coefficients(interior_coefficients), self.compute_points_from_coefficients(boundary_coefficients)
        points = i_points + b_points
        congruence_counter = 0
        for point in points:
            if self.square_of_size_of_vector(point) % D == j:
                congruence_counter += 1
                v1 = self.vi[0]
                v2 = self.vi[1]
                x_list = [point[0]+self.av1_plus_bv2(-0.5,0.5,v1,v2)[0],point[0]+self.av1_plus_bv2(-0.5,-0.5,v1,v2)[0],point[0]+self.av1_plus_bv2(0.5,-0.5,v1,v2)[0],point[0]+self.av1_plus_bv2(0.5,0.5,v1,v2)[0]]
                y_list = [point[1]+self.av1_plus_bv2(-0.5,0.5,v1,v2)[1],point[1]+self.av1_plus_bv2(-0.5,-0.5,v1,v2)[1],point[1]+self.av1_plus_bv2(0.5,-0.5,v1,v2)[1],point[1]+self.av1_plus_bv2(0.5,0.5,v1,v2)[1]]
                plt.fill(x_list, y_list, facecolor='lightsalmon', edgecolor='orangered')
            
        if tiled:
            tile_centres = self.lattice_points_contained_by_circle_of_radius_R(radius/D,False)
            tile_centres = [self.multiply_vector_by_scalar(tile_centre,D) for tile_centre in tile_centres]
            for point in tile_centres:
                s = (D)/2
                x_list = [point[0]-s,point[0]-s,point[0]+s,point[0]+s]
                y_list = [point[1]+s,point[1]-s,point[1]-s,point[1]+s]
                plt.fill(x_list, y_list, facecolor='none', edgecolor='green',linewidth=1.00)

        if (j == 0) and (D == 1):
            plt.scatter([point[0] for point in points],[point[1] for point in points],s=10/radius, color='lightseagreen') # maybe can make boundary points yellow to make visual nicer?
            plt.title(str(len(points))+' points', loc='left')
            plt.title('Radius ' + str(radius), loc='right')
            plt.show()
            return

        plt.scatter([point[0] for point in points],[point[1] for point in points],s=10/radius, color='red') # maybe can make boundary points yellow to make visual nicer?
        plt.title(str(len(points))+' points in circle, ' + str(congruence_counter) +' congruent to '+str(j) + ' mod '+str(D), loc='left')
        plt.title('Radius ' + str(radius), loc='right')
        plt.show()

    def is_prime(self,n):
        for prime in primes:
            if prime == n:
                return True
            if prime >= n:
                return False
        if n == 2 or n == 3:
            return True
        for i in range(2,int(sqrt(n))+2):
            if n % i == 0:
                return False
        return True

    def is_integer(self,n):
        return (n - floor(n) == 0)

    def is_gauss_prime(self,n):
        s = sqrt(n)
        if self.is_integer(s):
            s = int(s)
            return self.is_prime(s) and s % 4 == 3
        return self.is_prime(n) and n % 4 == 1

    def angle_of_point(self,point):
        x = point[0]
        y = point[1]
        if x == 0:
            if y > 0:
                return 0.5
            if y < 0:
                return -0.5
        else:
            v = atan(y/x)*1/pi
            if x < 0:
                if y < 0:
                    return v - 1
                if y >= 0:
                    return v + 1
            return v

    def gauss_primes_contained_in_sector_of_annulus(self,r1,r2,a,b):
        points = []
        for i in range(-r2,r2+1):
            valid_j = []
            for j in range(int(sqrt(r2**2-i**2))+1):
                if j != 0:
                    valid_j.append(j)
                    valid_j.append(-j)
                else:
                    valid_j.append(0)
            for j in valid_j:
                if r1**2 <= i**2 + j**2 <= r2**2:
                    if self.is_gauss_prime(i**2 + j**2):
                        if a <= self.angle_of_point([i,j]) <= b:
                            points.append([i,j])
        return points

    
    def plot_gauss_primes_contained_in_sector_of_annulus(self,r1,r2,a,b):
        """Plots gauss primes contained in a sector (angle a to b) of an annulus (radius r1 to r2)
        
        Parameters
        ----------
        r1 = annulus small radius
        r2 = annulus big radius
        a = first sector angle
        b = second sector angle

        Returns
        ------
        None
        """
        points = self.gauss_primes_contained_in_sector_of_annulus(r1,r2,a,b)
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

        # Major ticks every 10, minor ticks every 1
        major_ticks = np.arange(-2*r2 - (-2*r2)%10, 2*r2 -(2*r2)%10 +11, 10)
        minor_ticks = np.arange(-2*r2 - (-2*r2)%10, 2*r2 -(2*r2)%10 +11, 1)

        ax.set_xticks(major_ticks)
        ax.set_xticks(minor_ticks, minor=True)
        ax.set_yticks(major_ticks)
        ax.set_yticks(minor_ticks, minor=True)

        # And a corresponding grid
        ax.grid(which='both')

        # Or if you want different settings for the grids:
        ax.grid(which='minor', alpha=0.2)
        ax.grid(which='major', alpha=0.5)

        angle = [2*pi/1000*i for i in range(int(500*a),int(500*b)+1)]

        x = [r2*cos(i) for i in angle]
        y = [r2*sin(i) for i in angle]
        x2 = [r1*cos(i) for i in reversed(angle)]
        y2 = [r1*sin(i) for i in reversed(angle)]
        x = x2 + x + [x2[0]]
        y = y2 + y + [y2[0]]

        plt.axis("equal")
        plt.plot(x, y, color='black')
        points = self.compute_points_from_coefficients(points)
        for point in points:
            v1 = self.vi[0]
            v2 = self.vi[1]
            x_list = [point[0]+self.av1_plus_bv2(-0.5,0.5,v1,v2)[0],point[0]+self.av1_plus_bv2(-0.5,-0.5,v1,v2)[0],point[0]+self.av1_plus_bv2(0.5,-0.5,v1,v2)[0],point[0]+self.av1_plus_bv2(0.5,0.5,v1,v2)[0]]
            y_list = [point[1]+self.av1_plus_bv2(-0.5,0.5,v1,v2)[1],point[1]+self.av1_plus_bv2(-0.5,-0.5,v1,v2)[1],point[1]+self.av1_plus_bv2(0.5,-0.5,v1,v2)[1],point[1]+self.av1_plus_bv2(0.5,0.5,v1,v2)[1]]
            plt.fill(x_list, y_list, facecolor='lightsalmon', edgecolor='orangered')
        plt.scatter([point[0] for point in points],[point[1] for point in points],s=10/r2) # maybe can make boundary points yellow to make visual nicer?
        plt.title(str(len(points))+' points', loc='left')
        plt.title('Radii ' + str(r1)+','+str(r2), loc='right')
        plt.show()

f = open('primesuptomillionlist.txt', 'r')
for line in f.readlines():
    primes_string = line.split(', ')
global primes
primes = []
for prime_string in primes_string:
    primes.append(int(prime_string))
f.close()
latt = Lattice([1,0],[0,1])
latt.plot_lattice_points_contained_by_circle_of_radius_R(25)
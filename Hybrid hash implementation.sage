from sage.all import *
import time


###############################################################################
# Ask the user for a message and return its binary representation in a list
###############################################################################
def user_input():

    message = input("Enter your message: ")
    print("Typed message:", message)

    bits = [int(b) for c in message for b in format(ord(c), '08b')]
    return(bits)

###############################################################################
# Basic hash function H computation. Given curve points matrix P and a length 
# 'm' vector 'v', it computes the length 'n' product 'P*v' and returns a binary 
# string of length m' = n*(l+1) per result ('l' is the length of the x-coordinate 
# of elliptic curve points, 256 in the current implementation).
###############################################################################
def H(P,vi,n,m,PointInfinity,PointInfinitySequence,q):
    output=[]

    for i in range(n):

        # -- Compute the i-th matrix row - vector product
        acum = PointInfinity
        for j in range(m):
            acum = acum + vi[j]*P[i*m + j]

        # -- Take the resulting elliptic curve point x-coordinate or the established
        # -- representation of the point at infinity.
        
        if(acum == PointInfinity):
            bitlist = PointInfinitySequence
        else:
            bitlist = Integer(acum[0]).bits()

            # -- Prepend zeros if its length is below the maximum (bitlength of modulus 'q')
            while len(bitlist) < q.nbits():
                bitlist.insert(0,0)

            # -- Append the compressed y-coordinate
            if Integer(acum[1]) < (Integer(q)-1)/2 :
                bitlist.append(0)
            else:
                bitlist.append(1)

        output.extend(bitlist)

    return(output)

###############################################################################
# Full hash computation
###############################################################################
def HashComputation(P,vv,n,m,m_prime,PointInfinity,PointInfinitySequence,q):

    # -- Append zeros to the last block and compute its hash digest 
    numBlocks = len(vv)
    input = vv[numBlocks-1]
    while (len(input)<m):
        input.append(0)
    result = H(MatrixP,input,n,m,PointInfinity,PointInfinitySequence,q)

    # -- Include the rest of blocks in the computation of the hash digest 
    index = numBlocks-2
    while index >= 0:
        input=vv[index]
        input.extend(result)
        result = H(MatrixP,input,n,m,PointInfinity,PointInfinitySequence,q)
        index = index - 1
        
    return(result)


#####################################################
# Main procedure
#####################################################

# -- The lattice is defined by a matrix of 'n' x 'm' dimensions 
l = 256
n = 2
c = 1.5
m_prime = n*(l + 1)
m = ceil(c*m_prime)
d = m - m_prime

# -- Ask the user to type a string and return its binary representation in a list 
v = user_input()
print("Typed message is", len(v) ,"bits long.")

# -- Pad the binary string so that it can be divided into blocks of length 'd' 
# -- This padding includes a 64-bit sequence at the end representing the length 
# -- of the original message.

# -- Represent the input message length in a 64-bit sequence 
v_length = Integer(len(v)).bits()
while len(v_length) < 64 :
    v_length.append(0)

v_length.reverse()

# -- Pad the original message and append de 64-bit representation of original message length 
while len(v)%d != (d-64):
    v.append(0)
v.extend(v_length)
print("Padded message is", len(v) ,"bits long.")

# -- Split 'v' into 'v_1', 'v_2', ... , 'v_k' blocks 
vv=[]

n_blocks = (len(v)) / d
print("Number of blocks:",n_blocks)
print("Length of each block:",d,"bits\n")

for j in range(n_blocks):
     vj = v[j*d:(j+1)*d]
     vv.append(vj)
     print(f"v[{j+1}] =",vj)
print(" ")

# -- Elliptic curve setup 

# -- secp256k1: y^2 = x^3 + 7
q = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
F = GF(q)
E = EllipticCurve(F, [0, 7])
r = E.cardinality()
assert r.is_prime()
assert l == ceil(log(r, 2))
PointInfinity = E(0)

# -- We search for a value 'x' so that there does not exist a point '(x,y)' in curve 'E'.
# -- Such 'x' will be used to represent the point at infinity.
x=F(0)
while (x**3+F(7)).is_square():
    x = x + 1;
PointInfinitySequence = Integer(x).bits()
# -- Prepend zeros if its length is below that of 'q'. An extra zero is appended to get the length of 'q'+1 
while len(PointInfinitySequence) < q.nbits():
    PointInfinitySequence.insert(0,0)
PointInfinitySequence.append(0)

# -- Random generation of lattice matrix points 
start = time.time()
MatrixP = []
for i in range(n):
    for j in range(m):
        MatrixP.append(E.random_point())
print("Time to generate matrix of curve points:", time.time() - start, "seconds.\n")

# -- Hash digest computation
start = time.time()
result = HashComputation(MatrixP,vv,n,m,m_prime,PointInfinity,PointInfinitySequence,q)
print("Time to compute the hash digest:", time.time() - start, "seconds.")
print("Digest length:",len(result),"bits")
print(f"H*(v) =", result)

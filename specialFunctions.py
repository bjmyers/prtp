import numpy as np

def factorial(x):
    '''Inputs a list, returns a list of the factorial of each element of the input list
       Ignores negatives'''
    x = x.copy()
    ignore = (x < 0)
    y = np.where(ignore,0,x)
    y = np.array(list(map(np.math.factorial,y)))
    y = np.where(ignore,x,y)
    return y
    

def radialpoly(rho,n,m):
    output = 0
    if (rho<=1):
        for j in range((0,n-np.abs(n))/2):
            const = ((-1)**j)*np.math.factorial(n-j)/np.math.factorial(j)/np.math.factorial((n+m)/2-j)/factorial((n-m)/2-j)
            output = output + const*rho**(n-2*j)
    return output
    
    
def fastradpoly(rho,n,m):
    if(rho==0):
        if(m==0):
            return (-1.0)**(n/2)
        else:
            return 0
    
    if(m==n):
        return rho**n
    elif(m==n-2):
        return n*rho**n - (n-1)*rho**(n-2)
    else:
        tempm = n-4
        r4 = rho**n
        r2 = n*rho**n - (n-1)*rho**(n-2)
        while(tempm >= m):
            h3 = -4*(tempm+2)*(tempm+1)/(n+tempm+2)/(n-tempm)
            h2 = h3*(n+tempm+4)*(n-tempm-2)/4./(tempm+3) + (tempm+2)
            h1 = .5*(tempm+4)*(tempm+3) - (tempm+4)*h2 + h3*(n+tempm+6)*(n-tempm-4)/8.0
            
            if (tempm == m):
                return h1*r4 + (h2 + h3/rho**2)*r2
            
            r4 = r2
            r2 = h1*r4 + (h2 + h3/rho**2)*r2
            tempm -= 2


def fastradder(rho,n,m):
    if (rho==0):
        if(m==1):
            return (-1.0)**((n-1)/2)*(n+1)/2
        else:
            return 0
            
    if (m==n):
        return n*rho**(n-1)
    elif(m==n-2):
        return n*(n*rho**(n-1)) - (n-1)*(n-2)*rho**(n-3)
    else:
        tepmm = n-4
        r4 = n*rho**(n-1)
        r2 = n*(n*rho**(n-1)) - (n-1)*(n-2)*rho**(n-3)
        while(tempm >= m):
            h3 = -4*(tempm+2)*(tempm+1)/(n+tempm+2)/(n-tempm)
            h2 = h3*(n+tempm+4)*(n-tempm-2)/4.0/(tempm+3) + (tempm+2)
            h1 = .5*(tempm+4)*(tempm+3) - (tempm+4)*h2 + h3*(n+tempm+6)*(n-tempm-4)/8.0
            output = h1*r4 + (h2 + h3/rho**2)*r2 - 2*h3*fastradpoly(rho,n,tempm+2)/rho**3
            
            if(tempm == m):
                return output
            r4 = r2
            r2 = output
            tempm -= 2


def zernset(rho,theta,rorder,aorder,znum,polyout,derrho,dertheta):
    tznum,radnum = 1,1
    while(tznum < znum):
        tznum += radnum + 1
        radnum += 1
    rnm = np.zeros((radnum,radnum))
    rprime = np.zeros((radnum,radnum))
    
    for i in range(1,radnum):
        n = i-1
        for j in range(n+1,1,-2):
            m = j-1
            if (rho==0):
                if (m==1):
                    rprime[i,j] = (-1.0)**((n-1)/2)*(n+1)/2
                    rnm[i,j] = 0
                elif (m==0):
                    rprime[i,j] = 0
                    rnk[i,j] = (-1.0)**(n/2)
                else:
                    rprime[i,j] = 0
                    rnm[i,j] = 0
            else:
                if (n==m):
                    rprime[i,j] = n*rho**(n-1)
                    rnm[i,j] = rho**n
                elif (m==n-2):
                    rprime[i,j] = n*rprime[i,i] - (n-1)*rprime[i-2,i-2]
                    rnm[i,j] = n*rnm[i,i] - (n-1)*rnm[i-2,i-2]
                else:
                    h3 = -4*(m+2)*(m+1)/(n+m+2)/(n-m)
                    h2 = h3*(n+m+4)*(n-m-2)/4.0/(m+3) + (m+2)
                    h1 = .5*(m+4)*(m+3) - (m+4)*h2 + h3*(n+m+6)*(n-m-4)/8.0
                    rmm[i,j] = h1*rnm[i,j+4] + (h2+h3/rho**2)*rnm(i,j+2)
                    rprime[i,j] = h1*rprime[i,j+4] + (h2+h3/rho**2)*rprime[i,j+2] - 2*h3/rho**3*rnm[i,j+2]
                    
    for i in range(1,znum):
        n = rorder[i]
        mm = aorder[i]
        m = np.ans(mm)
        norm = np.sqrt(2*(n-1))
        if (mm<0):
            polyout[i] = norm * rnm[int(n)+1,int(m)+1] * np.sin(m*theta)
            derrho[i] = norm * rprime[int(n) + 1,int(m)+1] * np.sin(m*theta)
            dertheta[i] = norm * rnm[int(n)+1,int(m)+1] * np.cos(m*theta) * m
        elif (mm>0):
            polyout[i] = norm * rnm[int(n)+1,int(m)+1] * np.cos(m*theta)
            derrho[i] = norm * rprime[int(n)+1,int(m)+1] * np.cos(m*theta)
            dertheta[i] = -norm * rnm[int(n)+1,int(m)+1] * np.sin(m*theta) * m
        else:
            polyout[i] = norm*np.sqrt(0.5)*rnm[int(n)+1,int(m)+1]
            derrho[i] = norm*np.sqrt(0.5)*rprime[int(n)+1,int(m)+1]
            dertheta[i] = 0
    
    return polyout,derrho,dertheta


def zernike(rho,theta,n,m):
    rnm = fastradpoly(rho,float(n),np.abs(float(m)))
    
    norm = np.sqrt(2*(float(n)+1))
    if (m<0):
        return norm * rnm * np.sin(m*theta)
    elif (m>0):
        return norm * rnm * np.cos(m*theta)
    else:
        return norm*np.sqrt(0.5)*rnm


def radialder(rho,n,m):
    output = 0
    if(rho<=1):
        for j in range(0,(n-np.abs(m))/2):
            const = (-1)**j*np.math.factorial(n-j)/np.math.factorial(j)/np.math.factorial((n+m)/2-j)/np.math.factorial((n-m)/2-j)
            output += const*(n-2*j)*rho**(n-2*j-1)
    return output


def zernrhoder(rho,theta,n,m):
    rnm = fastradder(rho,float(n),np.abs(float(m)))
    norm = np.sqrt(2*(float(n)+1))
    if (m<0):
        return norm * rnm * np.sin(np.abs(m)*theta)
    elif (m>0):
        return norm * rnm * np.cos(m*theta)
    else:
        return norm*np.sqrt(0.5)*rnm


def zernthetader(rho,theta,n,m):
    rnm = fastradpoly(rho,float(n),np.abs(float(m)))
    norm = np.sqrt(2*(float(n)+1))
    if (m<0):
        return norm * rnm * np.abs(m) *  np.cos(np.abs(m)*theta)
    elif (m>0):
        return -1 * norm * rnm * m * np.sin(m*theta)
    else:
        return 0

def legendre(x,n):
    '''Takes two lists as inputs'''
    x = x.copy()
    n = n.copy()
    output = np.zeros(len(x))
    
    x = np.where(np.abs(x) > 1.0, x[:] / np.abs(x[:]),x[:])
    
    # Identify cases where no more calculations are needed
    done = (n == 0) 
    output = np.where(done,1,0)
    
    i = 0
    while True:
        needsaddition = ((np.floor(n[:])/2).astype(int) > i)
        alldone = np.logical_not(needsaddition).all()
        if (alldone):
            break
        output = np.where(np.logical_and(np.logical_not(done),needsaddition),output[:] + (-1)**(i)*factorial(2*n[:]-2*i)/np.math.factorial(i)/factorial(n[:]-i)/factorial(n[:]-2*i)/2**n[:]*x[:]**(n[:]-2*i),output[:])
        i += 1
        
    return output


def legendrep(x,n):
    '''Takes two lists as inputs'''
    x = x.copy()
    n = n.copy()
    output = np.zeros(len(x))
    
    case1 = (n==0)
    case2 = (n==1)
    output[case2] = 1
    case3 = np.logical_and((x==0),(n%2==0))
    
    i = 0
    while True:
        needsaddition = ((np.floor(n[:])/2).astype(int) > i)
        alldone = np.logical_not(needsaddition).all()
        if (alldone):
            break
        add1 = np.logical_and(np.logical_not(case1),np.logical_not(case2))
        add2 = np.logical_and(np.logical_not(case3),needsaddition)
        add = np.logical_and(add1,add2)
        output = np.where(add,output[:] + (-1)**(i)*factorial(2*n-2*i)/np.math.factorial(i)/factorial(n-i)/factorial(n-2*i)/2**n[:]*x[:]**(n-2*i),output[:])
        i += 1
    
    return output


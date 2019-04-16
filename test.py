from Rays import Rays
import numpy as np
import matplotlib.pyplot as plt

'''
test.py

This python script is used to test the functionality of Rays.py

Whenever prtp is used on a new machine or with a new version of Python, this
script should be run to make sure all functionality is retained
'''

# Keeps track of how many tests have passed
passed = 0
# Keeps track of how many tests there are total
total = 36
# Keeps track of messages that describe to the user which tests were failed
failedtests = []


## Rays Generation Tests

# Test 1 - Default Rays generation
try:
    x = np.array([1.3,1.2,1.1])
    y = np.array([2.1,0.0,3.0])
    z = np.array([.54,-.4,.89])
    l = np.array([14.,17.,.32])
    m = np.array([.1,3.45,1.2])
    n = np.array([7.8,5.6,.23])
    ux = np.array([-1.2,0,2.34])
    uy = np.array([10000,1e5,-.3])
    uz = np.array([1456,135.7,1.234])
    
    r = Rays(x,y,z,l,m,n,ux,uy,uz)
    if len(r) != 3:
        failedtests.append('Test 1 - default constructor did not successfully create an object of Length 3')
    else:
        passed += 1
except:
    failedtests.append('Test 1 - Default Rays construction method failure')


# Test 2 - pointsource Rays generation (0 angle)
try:
    r = Rays.pointsource(0,10)
    if len(r) != 10:
        failedtests.append('Test 2 - Rays.pointsource did not successfully create an object of Length 10')
    else:
        passed += 1
except:
    failedtests.append('Test 2 - Failure of Rays.pointsource when ang=0')


# Test 3 - pointsource Rays generation (Positive Angle)
try:
    r = Rays.pointsource(0.75,1)
    if len(r) != 1:
        failedtests.append('Test 3 - Rays.pointsource did not successfully create an object of Length 1')
    else:
        passed += 1
except:
    failedtests.append('Test 3 - Failure of Rays.pointsource when ang=.75')


# Test 4 - pointsource Rays generation (Negative Angle)
try:
    r = Rays.pointsource(-np.pi/2,100)
    if len(r) != 100:
        failedtests.append('Test 4 - Rays.pointsource did not successfully create an object of Length 100')
    else:
        passed += 1
except:
    failedtests.append('Test 4 - Failure of Rays.pointsource when ang=-np.pi/2')


# Test 5 - circularbeam Rays generation (0 radius)
try:
    r = Rays.circularbeam(0,103)
    if len(r) != 103:
        failedtests.append('Test 5 - Rays.circularbeam did not successfully create an object of Length 103')
    else:
        passed += 1
except:
    failedtests.append('Test 5 - Failure of Rays.circularbeam when rad=0')
    

# Test 6 - circularbeam Rays generation (Positive Radius)
try:
    r = Rays.circularbeam(3.5,25)
    if len(r) != 25:
        failedtests.append('Test 6 - Rays.circularbeam did not successfully create an object of Length 25')
    else:
        passed += 1
except:
    failedtests.append('Test 6 - Failure of Rays.circularbeam when rad=3.5')


# Test 7 - circularbeam Rays generation (Negative Radius)
try:
    r = Rays.circularbeam(-.36,3)
    if len(r) != 3:
        failedtests.append('Test 7 - Rays.circularbeam did not successfully create an object of Length 3')
    else:
        passed += 1
except:
    failedtests.append('Test 7 - Failure of Rays.circularbeam when rad=-.36')


# Test 8 - annulus Rays generation (Positive rout)
try:
    r = Rays.annulus(0,1.5,10)
    if len(r) != 10:
        failedtests.append('Test 8 - Rays.annulus did not successfully create an object of Length 10')
    else:
        passed += 1
except:
    failedtests.append('Test 8 - Failure of Rays.annulus when rin=0 and rout=1.5')


# Test 9 - annulus Rays generation (Negative rout)
try:
    r = Rays.annulus(0,-3,1)
    if len(r) != 1:
        failedtests.append('Test 9 - Rays.annulus did not successfully create an object of Length 1')
    else:
        passed += 1
except:
    failedtests.append('Test 9 - Failure of Rays.annulus when rin=0 and rout=-3')


# Test 10 - annulus Rays generation (specified zhat)
try:
    r = Rays.annulus(0,1.5,10,zhat=.5)
    if (r.n != .5).any():
        failedtests.append('Test 10 - Rays.annulus did not correctly specify zhat')
    elif len(r) != 10:
        failedtests.append('Test 10 - Rays.annulus did not successfully create an object of Length 10')
    else:
        passed += 1
except:
    failedtests.append('Test 10 - Failure of Rays.annulus when rin=0, rout=1.5, and zhat=.5')


# Test 11 - subannulus Rays generation (0 dphi)
try:
    r = Rays.subannulus(0,1.5,0,10)
    if len(r) != 10:
        failedtests.append('Test 11 - Rays.subannulus did not successfully create an object of Length 10')
    else:
        passed += 1
except:
    failedtests.append('Test 11 - Failure of Rays.subannulus when dphi=0')


# Test 12 - subannulus Rays generation (Positive dphi)
try:
    r = Rays.subannulus(0,-3,.5,1)
    if len(r) != 1:
        failedtests.append('Test 12 - Rays.subannulus did not successfully create an object of Length 1')
    else:
        passed += 1
except:
    failedtests.append('Test 12 - Failure of Rays.subannulus dphi=.5')


# Test 13 - subannulus Rays generation (specified zhat)
try:
    r = Rays.subannulus(0,1.5,-1.0,10,zhat=.5)
    if (r.n != .5).any():
        failedtests.append('Test 13 - Rays.subannulus did not correctly specify zhat')
    elif len(r) != 10:
        failedtests.append('Test 13 - Rays.subannulus did not successfully create an object of Length 10')
    else:
        passed += 1
except:
    failedtests.append('Test 13 - Failure of Rays.subannulus when dphi=-1.0 and zhat=.5')


# Test 14 - xslit Rays generation
try:
    r = Rays.xslit(0,-3,1,zhat=0.9)
    if len(r) != 1:
        failedtests.append('Test 14 - Rays.xslit did not successfully create an object of Length 1')
    else:
        passed += 1
except:
    failedtests.append('Test 14 - Failure of Rays.xslit')


# Test 15 - rectarray Rays generation
try:
    r = Rays.rectarray(12,10.5,10)
    if len(r) != 100:
        failedtests.append('Test 15 - Rays.rectarray did not successfully create an object of Length 100')
    else:
        passed += 1
except:
    failedtests.append('Test 15 - Failure of Rays.rectarray')


# Test 16 - convergingbeam Rays generation
try:
    r = Rays.convergingbeam(1,.5,1.5,0.0,0.5,25,.1)
    if len(r) != 25:
        failedtests.append('Test 16 - Rays.convergingbeam did not successfully create an object of Length 25')
    else:
        passed += 1
except:
    failedtests.append('Test 16 - Failure of Rays.convergingbeam')


# Test 17 - convergingbeam2 Rays generation
try:
    r = Rays.convergingbeam2(.5,0.1,0.3,-1.2,0.8,22,1.0)
    if len(r) != 22:
        failedtests.append('Test 17 - Rays.convergingbeam2 did not successfully create an object of Length 22')
    else:
        passed += 1
except:
    failedtests.append('Test 17 - Failure of Rays.convergingbeam2')


# Test 18 - rectbeam Rays generation
try:
    r = Rays.rectbeam(2,3,1)
    if len(r) != 1:
        failedtests.append('Test 18 - Rays.rectbeam did not successfully create an object of Length 1')
    else:
        passed += 1
except:
    failedtests.append('Test 18 - Failure of Rays.rectbeam')


## Tag Tests
x = np.array([1.3,1.2,1.1])
y = np.array([2.1,0.0,3.0])
z = np.array([.54,-.4,.89])
l = np.array([14.,17.,.32])
m = np.array([.1,3.45,1.2])
n = np.array([7.8,5.6,.23])
ux = np.array([-1.2,0,2.34])
uy = np.array([10000,1e5,-.3])
uz = np.array([1456,135.7,1.234])

r = Rays(x,y,z,l,m,n,ux,uy,uz)

# Test 19 - addTag
try:
    r.addTag('tag1',np.array([0,0,1]))
    passed += 1
except:
    failedtests.append('Test 19 - Rays.addTag failed')


# Test 20 - addTag (attempting a tagname starting with '~')
try:
    r.addTag('~tag2',np.array([1,2,3]))
    failedtests.append('Test 20 - Rays.addTag accepted a tag starting with "~"')
except ValueError:
    passed += 1


# Test 21 - addTag (attempting a tag with incorrect length)
try:
    r.addTag('tag2',np.array([1,2,3,4,5]))
    failedtests.append('Test 21 - Rays.addTag accepted a tag with incorrect length')
except Exception:
    passed += 1


# Test 22 - addTag (attempting an existing tagname)
try:
    r.addTag('tag1',np.array([1,2,3]))
    failedtests.append('Test 22 - Rays.addTag accepted a tag that already exists')
except Exception:
    passed += 1


# Test 23 - removeTag
try:
    r.removeTag('tag1')
    passed += 1
except:
    failedtests.append('Test 23 - Rays.removeTag failed')

# Test 24 - removeTag (removing a tag that doesn't exist)
try:
    r.removeTag('tag1')
    failedtests.append('Test 24 - Rays.removeTag removed a tag that doesnt exist')
except:
    passed += 1


# Test 25 - getTag
r.addTag('tag1',np.array([1,2,3]))
try:
    t = r.getTag('tag1')
    if (t==np.array([1,2,3])).all():
        passed += 1
    else:
        failedtests.append('Test 25 - Rays.getTag did not retrieve the correct Tag')
except:
    failedtests.append('Test 25 - Rays.getTag failed')


# Test 26 - getTag (getting a tag that doesn't exist)
try:
    t = r.getTag('tag99')
    if t is None:
        passed += 1
    else:
        failedtests.append('Test 26 - Rays.getTag retrieved a nonexistent Tag')
except:
    failedtests.append('Test 26 - Rays.getTag failed')


# Test 27 - combineTags
r.addTag('tag2',np.array([0,0,1]))
r.addTag('tag3',np.array([0,1,0]))
try:
    t = r.combineTags(['tag2','tag3'])
    if not (t==np.array([0,1,1])).all():
        raise Exception()
    t = r.combineTags('tag2 ~tag3',delim=" ")
    if not (t==np.array([1,0,1])).all():
        raise Exception()
    t = r.combineTags('tag2.tag3',delim=" ",orcombination=False)
    if not (t==np.array([0,0,0])).all():
        raise Exception()
    passed += 1
except:
    failedtests.append('Test 27 - Rays.combineTags failed')


## Param Tests
# Test 28 - addParam
try:
    r.addParam('Param1',np.array([0,0,1]))
    passed += 1
except:
    failedtests.append('Test 28 - Rays.addParam failed')


# Test 29 - addParam (attempting a Param with incorrect length)
try:
    r.addParam('Param2',np.array([1,2,3,4,5]))
    failedtests.append('Test 29 - Rays.addParam accepted a Param with incorrect length')
except Exception:
    passed += 1


# Test 30 - addParam (attempting an existing Paramname)
try:
    r.addParam('Param1',np.array([1,2,3]))
    failedtests.append('Test 30 - Rays.addParam accepted a Param that already exists')
except Exception:
    passed += 1


# Test 31 - removeParam
try:
    r.removeParam('Param1')
    passed += 1
except:
    failedtests.append('Test 31 - Rays.removeParam failed')

# Test 32 - removeParam (removing a Param that doesn't exist)
try:
    r.removeParam('Param1')
    failedtests.append('Test 32 - Rays.removeParam removed a Param that doesnt exist')
except:
    passed += 1


# Test 33 - getParam
r.addParam('Param1',np.array([1,2,3]))
try:
    t = r.getParam('Param1')
    if (t==np.array([1,2,3])).all():
        passed += 1
    else:
        failedtests.append('Test 33 - Rays.getParam did not retrieve the correct Param')
except:
    failedtests.append('Test 33 - Rays.getParam failed')


# Test 34 - getParam (getting a Param that doesn't exist)
try:
    t = r.getParam('Param99')
    if t is None:
        passed += 1
    else:
        failedtests.append('Test 34 - Rays.getParam retrieved a nonexistent Param')
except:
    failedtests.append('Test 34 - Rays.getParam failed')



## Addition Tests
x = np.array([1.3,1.2,1.1]) * -1
y = np.array([2.1,0.0,3.0]) * -1
z = np.array([.54,-.4,.89]) * -1
l = np.array([14.,17.,.32]) * -1
m = np.array([.1,3.45,1.2]) * -1
n = np.array([7.8,5.6,.23]) * -1
ux = np.array([-1.2,0,2.34]) * -1
uy = np.array([10000,1e5,-.3]) * -1
uz = np.array([1456,135.7,1.234]) * -1

r2 = Rays(x,y,z,l,m,n,ux,uy,uz)

# Test 35 - Basic Addition
try:
    r3 = r + r2
    if not (r3.x==np.array([1.3,1.2,1.1,-1.3,-1.2,-1.1])).all():
        failedtests.append('Test 35 - Rays + Rays failed to create correct x')
        raise Exception()
    if not (r3.y==np.array([2.1,0.0,3.0,-2.1,-0.0,-3.0])).all():
        failedtests.append('Test 35 - Rays + Rays failed to create correct y')
        raise Exception()
    if not (r3.z==np.array([.54,-.4,.89,-.54,.4,-.89])).all():
        failedtests.append('Test 35 - Rays + Rays failed to create correct z')
        raise Exception()
    if not (r3.l==np.array([14.,17.,.32,-14.,-17.,-.32])).all():
        failedtests.append('Test 35 - Rays + Rays failed to create correct l')
        raise Exception()
    if not (r3.m==np.array([.1,3.45,1.2,-.1,-3.45,-1.2])).all():
        failedtests.append('Test 35 - Rays + Rays failed to create correct m')
        raise Exception()
    if not (r3.n==np.array([7.8,5.6,.23,-7.8,-5.6,-.23])).all():
        failedtests.append('Test 35 - Rays + Rays failed to create correct n')
        raise Exception()
    if not (r3.ux==np.array([-1.2,0,2.34,1.2,-0,-2.34])).all():
        failedtests.append('Test 35 - Rays + Rays failed to create correct ux')
        raise Exception()
    if not (r3.uy==np.array([10000,1e5,-.3,-10000,-1e5,.3])).all():
        failedtests.append('Test 35 - Rays + Rays failed to create correct uy')
        raise Exception()
    if not (r3.uz==np.array([1456,135.7,1.234,-1456,-135.7,-1.234])).all():
        failedtests.append('Test 35 - Rays + Rays failed to create correct uz')
        raise Exception()
    if not (r3.getTag('tag1')==np.array([1,2,3,0,0,0])).all():
        failedtests.append('Test 35 - Rays + Rays failed to create correct tag1')
        raise Exception()
    if not (r3.getParam('Param1')==np.array([1,2,3,0,0,0])).all():
        failedtests.append('Test 35 - Rays + Rays failed to create correct Param1')
        raise Exception()
    passed += 1
except:
    failedtests.append('Test 35 - Rays + Rays failed')


# Test 36 - Augmented Addition Assignment
try:
    r2 += r
    if not (r2.x==np.array([-1.3,-1.2,-1.1,1.3,1.2,1.1])).all():
        failedtests.append('Test 36 - Rays += Rays failed to create correct x')
        raise Exception()
    if not (r2.y==np.array([-2.1,-0.0,-3.0,2.1,0.0,3.0])).all():
        failedtests.append('Test 36 - Rays += Rays failed to create correct y')
        raise Exception()
    if not (r2.z==np.array([-.54,.4,-.89,.54,-.4,.89])).all():
        failedtests.append('Test 36 - Rays += Rays failed to create correct z')
        raise Exception()
    if not (r2.l==np.array([-14.,-17.,-.32,14.,17.,.32])).all():
        failedtests.append('Test 36 - Rays += Rays failed to create correct l')
        raise Exception()
    if not (r2.m==np.array([-.1,-3.45,-1.2,.1,3.45,1.2])).all():
        failedtests.append('Test 36 - Rays += Rays failed to create correct m')
        raise Exception()
    if not (r2.n==np.array([-7.8,-5.6,-.23,7.8,5.6,.23])).all():
        failedtests.append('Test 36 - Rays += Rays failed to create correct n')
        raise Exception()
    if not (r2.ux==np.array([1.2,-0,-2.34,-1.2,0,2.34])).all():
        failedtests.append('Test 36 - Rays += Rays failed to create correct ux')
        raise Exception()
    if not (r2.uy==np.array([-10000,-1e5,.3,10000,1e5,-.3])).all():
        failedtests.append('Test 36 - Rays += Rays failed to create correct uy')
        raise Exception()
    if not (r2.uz==np.array([-1456,-135.7,-1.234,1456,135.7,1.234])).all():
        failedtests.append('Test 36 - Rays += Rays failed to create correct uz')
        raise Exception()
    if not (r2.getTag('tag1')==np.array([0,0,0,1,2,3])).all():
        failedtests.append('Test 36 - Rays += Rays failed to create correct tag1')
        raise Exception()
    if not (r2.getParam('Param1')==np.array([0,0,0,1,2,3])).all():
        failedtests.append('Test 36 - Rays += Rays failed to create correct Param1')
        raise Exception()
    passed += 1
except:
    failedtests.append('Test 36 - Rays += Rays failed')
    
















print('Passed {:d}/{:d} ({:.3f}%) tests'.format(passed,total,(passed/total*100)))
if len(failedtests) > 0:
    print('Failed Tests:')
    for f in failedtests:
        print(f)














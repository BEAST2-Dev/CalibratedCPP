import math

def CDF(birth_rate, death_rate, rho, time):
    num = rho * birth_rate * (1 - math.exp(-time * (birth_rate - death_rate)))
    denom = rho * birth_rate + (birth_rate * (1 - rho) - death_rate) * math.exp(-time * (birth_rate - death_rate))
    Q = num / denom
    return Q

def cherryProb(birth_rate, death_rate, rho, time, cherryTime):
    CDFt = CDF(birth_rate=birth_rate, death_rate=death_rate, rho=rho, time=time)
    CDFx = CDF(birth_rate=birth_rate, death_rate=death_rate, rho=rho, time=cherryTime)
    probLoc = (CDFt-CDFx)/(2*(CDFt-CDFx)+(CDFt-CDFx)**2)
    prob = probLoc * (CDFt + CDFx)/(2*CDFt) * 2
    return prob

def cherryProb_root(birth_rate, death_rate, rho, time, cherryTime):
    CDFt = CDF(birth_rate=birth_rate, death_rate=death_rate, rho=rho, time=time)
    CDFx = CDF(birth_rate=birth_rate, death_rate=death_rate, rho=rho, time=cherryTime)
    probLoc = 1/(1+2*(CDFt-CDFx))
    return probLoc


#fxn to pull simulation params out of file names
def parse_params(fn):
    x=fn.split('_')
    opt=x[3]
    mu=x[5]
    return {'opt':opt,'mu':mu,'sigmu':0.25}

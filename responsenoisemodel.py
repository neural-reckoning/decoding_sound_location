from shared import *
from numpy.random import poisson

class ResponseNoiseModel(object):
    '''
    Base class for different types of response noise. The __call__ method
    should return the results of applying the response noise to the response.
    The user can also specify a function rather than a class.
    '''
    def __call__(self, response):
        return response
    def __str__(self):
        return 'NoResponseNoise'
    __repr__ = __str__


NoResponseNoise = ResponseNoiseModel


class IndependentWhiteResponseNoise(ResponseNoiseModel):
    def __init__(self, sigma, rectify=True):
        self.sigma = sigma
        self.rectify = rectify
        
    def __call__(self, response):
        response = response+randn(len(response))*self.sigma
        if self.rectify:
            response[response<0] = 0
        return response      

    def __str__(self):
        return 'IndependentWhiteResponseNoise(%g,%s)'%(self.sigma, str(self.rectify))
    __repr__ = __str__


class PoissonResponseNoise(ResponseNoiseModel):
    def __call__(self, response):
        response = poisson(response)
        return response

    def __str__(self):
        return 'PoissonResponseNoise()'
    __repr__ = __str__

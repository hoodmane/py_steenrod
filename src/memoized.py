import collections
import functools

class memoized(object):
   """Decorator. Caches a function's return value each time it is called.
   If called later with the same arguments, the cached value is returned
   (not reevaluated).
   """
   def __init__(self, func):
      self.func = func
      self.cache = {}
      functools.update_wrapper(self, func)
      
   def __call__(self, *args, **kw_args):
      key = args + tuple(kw_args.items())
      try:
         key in self.cache
      except TypeError:
         return self.func(*args, **kw_args)
      if key in self.cache:
         #print("hit ", key)
         return self.cache[key]
      else:
         #print("miss ", key)
         value = self.func(*args, **kw_args)
         self.cache[key] = value
         return value

   def __get__(self, obj, objtype):
      '''Support instance methods.'''
      return functools.partial(self.__call__, obj)

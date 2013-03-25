#!/usr/bin/env python
# simple optimization algorithm for tuning parameters using UCT
from random import uniform, randrange, sample
import subprocess
from math import sqrt, log

APP_PATH="./test_uct_optimizer.py"

RUN_COMMAND=APP_PATH

class IntervalNode(object):
	"""
	Simple node representing an interval that we are exploring in the tree.
	TODO: use properties instead of so many getters.
	The intervals are not relative to their parents but instead are absolute numerical
	values of the original interval. As the tree becomes deeper we split the intervals into 
	finer and finer sub intervals.
	"""
	WIN_RESULTS = {'W': 1, 'L': 0, 'D': 0}
	

	def __init__(self, interval, uct_coeff, num_chunks):
		"""
		Simple tree node that allows us to iteratively refine
		the interval.
		The interval chunks sort of correspond to actions and a particular
		interval corresponds to states.
		"""
		self.children 		= [None for i in range(num_chunks)] # initialize to none
		self.interval	    = interval
		self.wins 			= 0
		self.visits 		= 0
		self.samples 		= []
		self.parent 		= None
		self.uct_coeff 		= uct_coeff
		self.num_chunks 	= num_chunks
		
	def randomSubinterval(self):
		"""
		Returns a random sub-interval as a tuple. The first parameter
		is the index of the left hand side of the chunk
		"""
		
		# get the left hand side starting range
		chunk = randrange(self.num_chunks)
		
		return (chunk, chunk * self.chunk_width, (chunk + 1) * self.chunk_width)
		
	def isRoot(self):
		"""Checks whether or not we are a root node by checking if we have an interval set"""
		return self.interval is None or self.parent is None
	
	def inflateRandomChild(self):
		"""
		Should be called when we have not examined all the intervals.
		One option might be to inflate all the children but only set
		their visits count to zero?
		"""
		# if we have inflated all the children don't continue
		if self.visits >= self.num_chunks:
			return None
			
		print "Visits: %d" % (self.visits, )
		
		# create a list of indicies of unvisited nodes
		unvisited = []
		
		for index in range(self.num_chunks):
			if self.children[index] is None:
				unvisited.append(index)
				
		# pick a random unvisisted sub_interval index
		rand_chunk_index = sample(unvisited, 1)[0]
		
		#if a child doesn't already exist, we make a new one
		node = self.makeChildForInterval(rand_chunk_index)
			
		return node
		
	def bestChild(self):
		"""Finds the best child representing an actino"""
		if len(self.children) == 0:
			return None
			
		if self.visits < self.num_chunks:
			# we should hit each kid at least once.
			# grab one of the unvisited children.
			# the calling code should then call inflateNewChild
			return None
		else:
			print repr([child.parent for child in self.children])
			# if we have visited all the children, then we should simply take the one with the highest UCT score
			return max(self.children, key=lambda child: child.uctValue)
	
	def findBestInterval(self):
		"""Descends the tree picking the best child at every step"""
		cur_node = self
		
		# just follow the best nodes down the root
		while cur_node.hasChildren:
			cur_node = cur_node.bestChild()
			
		return cur_node.interval
	
	@property
	def interval_width(self):
		"""Gets the total interval width"""
		return self.interval[1] - self.interval[0]
		
	@property
	def chunk_width(self):
		"""Gets the width of each subinterval"""
		return self.interval_width / self.num_chunks
		
	def makeChildForInterval(self, interval_index):
		"""Adds a new child if it doesn't already exist for the given interval index"""
		if self.children[interval_index] is not None:
			return self.children[interval_index]
			
		self.children[interval_index] = IntervalNode((interval_index * self.chunk_width, (interval_index + 1) * self.chunk_width), self.uct_coeff, self.num_chunks)
		
		# now add the child to ourselves
		self.addChild(self.children[interval_index], interval_index)
		
		return self.children[interval_index]
		
	def hasChildForInterval(self, interval):
		"""Checks to see if we have inflated a child for a given interval index"""
		return self.children[interval] is None
	
	@property	
	def uctValue(self):
		""" Computes the UCT ranking for this node. Should only be called on children node (never the root)"""
		if self.isRoot:
			return None
			
		# compute the UCT value. We bias our win rate based on our parents visits and our own visits
		return (self.wins / self.visits) + self.uct_coeff * sqrt(log(self.parent.visits) / self.visits)
	
	@property
	def hasChildren(self):
		return len(self.children) > 0
	
	@property
	def isLeaf(self):
		return not self.hasChildren
		
	def __repr__(self):
		return '[%.5f, %.5f]' % (self.interval[0], self.interval[1])
		
	def propegateResult(self, result):
		"""Propegates the UCT results back up the tree. The tuple should be the parameter, the value, and the result."""
		self.visits += 1
		
		# update based on result
		self.wins += self.WIN_RESULTS[result[2]]
		self.samples.append(result)
		
		# propegates the result up the tree
		if self.parent is not None:
			self.parent.propegateResult(result)
			
	def addChild(self, child, chunk_index):
		"""Sets the child at the specified chunk index"""
		self.children[chunk_index] = child
		child.parent = self
			
class UCTOptimizer(object):
	"""Uses UCT to optimize the values of a number of parameters"""
	
	def __init__(self, test_program, iterations, uct_coeff):
		self.test_program = test_program
		self.iterations   = iterations
		self.uct_coeff	  = uct_coeff
		self.num_chunks   = 10 # number of chunks on each level
		self.horizen 	  = 20 # the depth of the tree. Ideally we should calculate this based on the desired precision
		
	def sample(self, node, param):
		"""Takes a 'sample' by running a game using the parameters within this node's interval"""
		
		# take a uniform point from the interval
		p = uniform(node.interval[0], node.interval[1])
		
		# concatenate the list of parameters with keys and values
		param_str = '%s %.5f' % (param, p)
		
		# build the command for sampling
		command = self.test_program + ' ' + param_str
		
		try:
			win_loss = subprocess.check_output(command.split(' ')).strip()
		except subprocess.CalledProcessError as e:
			win_loss = e.output.strip()
			
		print 'Sampled: %s' % (win_loss, )
		# inform the node about the samples
		node.propegateResult((param, p, win_loss))
		
		return win_loss
		 
	def tuneParams(self, params):
		
		"""
		Tunes parameters should be list of tuples (name, start, end) against the program.
		Returns the tuned parameters.
		"""
		# start at the root and just work our way down
		# split into top level interval nodes
		# right now we only use one parameter but we'd like to use others in the future
		p1 		 = params[0]
		
		chunk = (p1[2] - p1[1]) / self.num_chunks
		
		root = IntervalNode((p1[1], p1[2]), self.uct_coeff, self.num_chunks)
		
		# now we simply start exploring the tree
		# do we pick a random arm or do we explore the arm further?
		cur_node 	= root
		parent_node = root
		
		# run up and down the tree the specified number of iterations
		for iterations in xrange(self.iterations):
			cur_node = parent_node = root
			
			print "Iteration %d" % (iterations, )
			# run down the tree until we hit a node without all of it's children expanded
			# or we hit a leaf node and want to expand the tree
			while cur_node is not None:
				
				print "heading down the tree"
				# the best child call will also tend to inflate other children at the same level
				# it returns None if we need to inflate a child
				parent_node = cur_node
				cur_node = parent_node.bestChild()
				
				
			# oops! we've fallen out of the tree. Randomly add a new if we are at a leaf node with a new sample.
			# no notion of "playouts" unless we thought about fixed tree depth
			new_child = parent_node.inflateRandomChild()
			
			self.sample(new_child, p1[0])
		
		# now find the best interval
		print root.findBestInterval()
		
if __name__ == "__main__":
	tuner = UCTOptimizer(RUN_COMMAND, 100, .3)
	params = tuner.tuneParams([('k', 0.0, 1000.0)])
	
	print params
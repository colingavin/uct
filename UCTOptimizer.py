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
		
	def random_subinterval(self):
		"""
		Returns a random sub-interval as a tuple. The first parameter
		is the index of the left hand side of the chunk
		"""
		
		# get the left hand side starting range
		chunk = randrange(self.num_chunks)
		
		return (chunk, chunk * self.chunk_width, (chunk + 1) * self.chunk_width)
	
	def is_root(self):
		"""Checks whether or not we are a root node by checking if we have an interval set"""
		return self.interval is None or self.parent is None
	
	def inflate_random_child(self):
		"""
		Should be called when we have not examined all the intervals.
		One option might be to inflate all the children but only set
		their visits count to zero?
		"""
		# if we have inflated all the children don't continue
		if self.visits > self.num_chunks:
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
		node = self.make_child_for_interval(rand_chunk_index)
			
		return node
		
	def best_child(self):
		"""Finds the best child representing an actino"""
		if len(self.children) == 0:
			return None
		
		if len(filter(lambda child: child is None, self.children)) > 0:
			# we should hit each kid at least once.
			# grab one of the unvisited children.
			# the calling code should then call inflateNewChild
			return None
		else:
			# if we have visited all the children, then we should simply take the one with the highest UCT score
			return max(self.children, key=lambda child: child.uct_value)
	
	def most_visited_child(self):
		children = filter(lambda child: child is not None, self.children)
		if len(children) == 0:
			return None
		return max(children, key=lambda child: child.visits)

	def find_best_interval(self):
		"""Descends the tree picking the best child at every step"""
		cur_node = self
		
		# just follow the best nodes down the root
		while cur_node.has_children:
			best_child = cur_node.most_visited_child()
			if best_child is not None:
				cur_node = best_child
			else:
				break

		return cur_node.interval
	
	@property
	def interval_width(self):
		"""Gets the total interval width"""
		return self.interval[1] - self.interval[0]
		
	@property
	def chunk_width(self):
		"""Gets the width of each subinterval"""
		return self.interval_width / self.num_chunks
		
	def make_child_for_interval(self, interval_index):
		"""Adds a new child if it doesn't already exist for the given interval index"""
		if self.children[interval_index] is not None:
			return self.children[interval_index]
			
		self.children[interval_index] = IntervalNode((self.interval[0] + interval_index * self.chunk_width, self.interval[0] + (interval_index + 1) * self.chunk_width), self.uct_coeff, self.num_chunks)
		
		# now add the child to ourselves
		self.add_child(self.children[interval_index], interval_index)
		
		return self.children[interval_index]
		
	def has_child_for_interval(self, interval):
		"""Checks to see if we have inflated a child for a given interval index"""
		return self.children[interval] is None
	
	@property
	def uct_value(self):
		""" Computes the UCT ranking for this node. Should only be called on children node (never the root)"""
		if self.is_root():
			return None
			
		# compute the UCT value. We bias our win rate based on our parents visits and our own visits
		return (float(self.wins) / float(self.visits)) + self.uct_coeff * sqrt(log(self.parent.visits) / self.visits)
	
	@property
	def has_children(self):
		return len(self.children) > 0
	
	@property
	def is_leaf(self):
		return not self.has_children

	def __repr__(self):
		def repr_indented(node, tabs):
			children = ""
			for child in node.children:
				if child is not None:
					children += tabs + repr_indented(child, tabs + '\t') + '\n'
			return '%i/%i@<%.5f, %.5f> :\n%s' % (node.wins, node.visits, node.interval[0], node.interval[1], children)
		return repr_indented(self, '\t')
		
	def propegate_result(self, result):
		"""Propegates the UCT results back up the tree. The tuple should be the parameter, the value, and the result."""
		self.visits += 1
		
		# update based on result
		self.wins += self.WIN_RESULTS[result[2]]
		self.samples.append(result)
		
		# propegates the result up the tree
		if self.parent is not None:
			self.parent.propegate_result(result)
			
	def add_child(self, child, chunk_index):
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
		self.samples_loss = []
		self.samples_win  = []
		
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
		
		if win_loss == 'W':
			self.samples_win.append(p)
		else:
			self.samples_loss.append(p)

		print 'Sampled: %s' % (win_loss, )
		# inform the node about the samples
		node.propegate_result((param, p, win_loss))
		
		return win_loss

	def output_mathematica(self, l):
		c = "{"
		for idx, sample in enumerate(l):
			c += str(sample)
			if idx != len(l) - 1:
				c += ", "
		c += "};"
		return c

	def tune_params(self, params):
		
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
				cur_node = parent_node.best_child()

			# oops! we've fallen out of the tree. Randomly add a new if we are at a leaf node with a new sample.
			# no notion of "playouts" unless we thought about fixed tree depth
			new_child = parent_node.inflate_random_child()
			
			self.sample(new_child, p1[0])
		
		print root

		print self.output_mathematica(self.samples_win)
		print self.output_mathematica(self.samples_loss)

		# now find the best interval
		best_interval = root.find_best_interval()
		return (params[0][0], best_interval[0], best_interval[1])
		
if __name__ == "__main__":
	tuner = UCTOptimizer(RUN_COMMAND, 400, 0.3)
	params = tuner.tune_params([('k', 0.0, 1000.0)])
	
	print params

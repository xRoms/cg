from sortedcontainers import SortedListWithKey
from test_parameters import bound, number_of_segments, seed
from random import uniform
from structs import rounded, Point, Pixel, Segment, SweepLine as SL
import numpy as np

starting_segments = []

#starting_segments = [Segment(Point(0, 3, 1, homogeneous=True), Point(48, 30, 10, homogeneous=True)),
#					 Segment(Point(49, 26, 10, homogeneous=True), Point(6, 8, 1, homogeneous=True)),
#					 Segment(Point(3, 8, 1, homogeneous=True), Point(6, 2, 1, homogeneous=True))]

hot = []

def generate_segs(n, seed = None):
	if (seed is not None):
		np.random.seed(seed)

	rand = lambda : np.random.randint(0, 5 * bound)
	for i in range(n):
		q, w= rand(), rand()
		while q == w:
			w = rand()
		starting_segments.append(Segment(Point(q, rand(), 5, homogeneous=True), Point(w, rand(), 5, homogeneous=True)))

generate_segs(number_of_segments, seed)

for q in starting_segments:
	print(q.start, " ", q.end)

current = SortedListWithKey(key=lambda pix: pix.center.y)

line = SL(starting_segments)

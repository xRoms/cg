from __future__ import print_function
import sys
import numpy as np
from enum import Enum
from sortedcontainers import SortedListWithKey, SortedDict
import heapq
from numpy.linalg import det
from test_parameters import bound, eps, handlers
from cg import Point, turn
from cg.utils import gcd
from collections import deque
from blist import blist

pixels = {}

segsinstatus = {}

pixelspassed = {}

def average(ceps, point):
	"""
	приводит данную точность и точку на "один уровень" однородных координат 

	:param ceps: float
	:param point: экземпляр Point
	:return: (int, Point)
	"""

	zeps = 1
	while (((zeps * ceps)%1 != 0)):
		zeps *= 10
	m = gcd(point.z, zeps)
	new_eps = int(zeps * ceps * (point.z // m))
	new_point = Point(int(point.x * (zeps // m)), int(point.y * (zeps // m)), int(point.z * (zeps // m)), homogeneous=True)
	return (new_eps, new_point)

def point_inside(pixel, segment):
	"""
	возвращает точку, принадлежащую segment, находящуюся внутри pixel

	:param pixel: экземпляр Pixel 
	:param segment: экземпляр Segment
	:return: экземпляр Point
	"""

	tmp = []
	cnt = 0
	for q in segment.intersections(pixel):
		if (q is not None):
			tmp.append(q)
			cnt += 1
	if cnt == 0:
		return None
	if cnt >= 2:
		temp = tmp[0] + tmp[1]	
		return Point(int(temp.x), int(temp.y), int(temp.z * 2), homogeneous=True)
	return tmp[0]

def add_pixel_to_seg(pixel, segment):
	"""
	добавляет запись, что segment пересекает pixel в pixelpassed

	:param pixel: экземпляр Pixel
	:param segment: экземпляр Segment
	"""

	if (segment in pixelspassed) and (pixelspassed[segment] != None):
		pixelspassed[segment].append((point_inside(pixel, segment), pixel.center))
	else:
		pixelspassed[segment] = [(point_inside(pixel, segment), pixel.center)]


def get_pixel(a):
	"""
	возвращает экземпляр Pixel по точке в нем

	:param a: экземпляр Point
	:return: экземпляр Pixel
	"""

	na = rounded(a)
	if (na not in pixels):
		pixels[na] = Pixel(na)
	return pixels[na]

def smart_gcd(a, b):
	"""
	наименьший общий делитель, игнорирующий нули

	:param a, b: int
	:return: int
	"""

	if a == 0:
		return b
	if b == 0:
		return a
	return gcd(a, b)

def normalize(smth):
	"""
	возвращает точку с упрощенными однородными координатами

	:param smth: экземпляр Point
	:return: экземпляр Point
	"""

	m = smart_gcd(smth.x, smart_gcd(smth.y, smth.z))
	q = 1
	if (smth.z <= 0):
		q = 1
	return Point(int(q * (smth.x // m)), int(q * (smth.y // m)), int(q * (smth.z // m)), homogeneous=True)

class Segment:
	"""
	Класс отрезка
	"""

	def __init__(self, a, b, isbound=0):
		"""
		Конструктор отрезка

		:param a, b: концы отрезка, экземпляры Point
		:param isbound: является ли отрезок верхней (-1), нижней (1) границей или же нет (0)
		"""

		self.start = normalize(min(a, b))
		self.end = normalize(max(a, b))
		self.isbound = isbound

	def __eq__(self, other):
		"""
		проверяет два отрезка на равенство

		:param other: экземпляр Segment
		:return: bool
		"""

		return self.start == other.start and self.end == other.end		

	def __hash__(self):

		"""
		хэш отрезка
		"""

		return hash((self.start.x, self.start.y, self.start.z, self.end.x, self.end.y, self.end.z))

	def atX(self, x):
		"""
		вычисляет y координату отрезка при данном x

		:param x: экземпляр Segment
		:return: числовое значение
		"""

		sx = self.start.x / self.start.z
		ex = self.end.x / self.end.z - sx
		dx = x - sx
		sy = self.start.y / self.start.z
		ey = self.end.y / self.end.z - sy
		sy += self.isbound * 0.00001
		if (dx > ex):
			return round(sy, 6)
		return round(dx/ex * ey + sy, 6)

	def intersects(self, other):
		"""
		пересечение отрезка с другим отрезком

		:param other: экземпляр Segment
		:return: Point точка пересечения или None, если пересечения нет
		"""

		a, b = normalize(self.start), normalize(self.end)
		c, d = normalize(other.start), normalize(other.end)
		if (a == c) and (b == d):
			return None
		acd, bcd = turn(a, c, d), turn(b, c, d)
		cab, dab = turn(c, a, b), turn(d, a, b)

		do_intersect = False
		if acd == bcd == cab == dab == 0:
			do_intersect = (a <= c <= b or a <= d <= b or c <= a <= d or c <= b <= d)
		else:
			do_intersect = acd != bcd and cab != dab
		
		if do_intersect:
			cross = lambda a, b: np.cross([a.coord], [b.coord])
			prod = np.cross(cross(a, b), cross(c, d))[0]
			if (prod[0] == 0 and prod[1] == 0 and prod[2] == 0):
				return max(b,d)
			return normalize(Point(prod, homogeneous = True))
		return None

	def intersections(self, pixel):
		"""
		возвращает список пересечений отрезка с каждой стороной пикселя

		:param pixel: экземпляр Pixel
		"""

		return [self.intersects(pixel.right()), self.intersects(pixel.bottom()), self.intersects(pixel.top()), self.intersects(pixel.left())]

	def __str__(self):
		return '({}, {})'.format(self.start, self.end)

def rounded(point):
	"""
	округляет точку до центра ближайшего пикселя

	:param x: экземпляр Point
	:return: Point, являющийся центром пикселя, которому принадлежит x
	"""

	res = average(eps, point)
	new_eps = res[0]
	new_point = res[1]
	return normalize(Point(int(halfround(new_point.x, new_eps)), int(halfround(new_point.y, new_eps)), int(new_point.z), homogeneous= True))


def halfround(a, new_eps):
	"""
	округляет число с определенной точностью

	:param a: числовое значение
	:param new_eps: числовое значение
	:return: числовое значение 
	"""
	if (a % new_eps <= new_eps / 2): 
		return a - a % new_eps
	else:
		return a + new_eps - a % new_eps 


class Pixel:
	"""
	Класс пикселя
	"""

	def __init__(self, point):
		"""
		Конструктор пикселя

		:param point: экземпляр Point, точка внутри пикселя
		"""

		self.center = rounded(point)
		self.lower = []
		self.upper = []
		self.segs = []

	def __eq__(self, other):
		"""
		проверяет два пикселя на равенство

		:param other: экземпляр Pixel
		:return: bool
		"""

		return self.center == other.center

	def __lt__(self, other):
		"""
		проверяет меньше ли один пиксель другого

		:param other: экземпляр Pixel
		:return: bool
		"""

		return self.center < other.center

	

	def top(self):
		"""
		возвращает отрезок, соответствующий верхней грани пикселя

		:return: экземпляр Segment
		"""

		return Segment(self.nw, self.ne, isbound = -1)
	
	def bottom(self):
		"""
		возвращает отрезок, соответствующий нижней грани пикселя

		:return: экземпляр Segment
		"""

		return Segment(self.sw, self.se, isbound = 1)

	def left(self):
		"""
		возвращает отрезок, соответствующий левой грани пикселя

		:return: экземпляр Segment
		"""

		return Segment(self.nw, self.sw)

	def right(self):
		"""
		возвращает отрезок, соответствующий правой грани пикселя

		:return: экземпляр Segment
		"""

		return Segment(self.ne, self.se)

	def is_on_top(self, point):
		"""
		проверяет, лежит ли точка на верхней грани пикселя

		:param point: экземпляр Point
		:return: bool
		"""

		new_self, new_point = self.center.same_level(point)
		yself = new_self[1]
		ypoint = new_point[1]
		cureps = eps / 2 * new_self[2]

		return ypoint == yself + cureps

	def is_on_bottom(self, point):
		"""
		проверяет, лежит ли точка на нижней грани пикселя

		:param point: экземпляр Point
		:return: bool
		"""

		new_self, new_point = self.center.same_level(point)
		yself = new_self[1]
		ypoint = new_point[1]
		cureps = eps / 2 * new_self[2]

		return ypoint == yself - cureps


	def is_on_left(self, point):
		"""
		проверяет, лежит ли точка на левой грани пикселя

		:param point: экземпляр Point
		:return: bool
		"""

		new_self, new_point = self.center.same_level(point)
		xself = new_self[0]
		xpoint = new_point[0]
		cureps = eps / 2 * new_self[2]

		return xpoint == xself - cureps

	def is_on_right(self, point):
		"""
		проверяет, лежит ли точка на правой грани пикселя

		:param point: экземпляр Point
		:return: bool
		"""

		new_self, new_point = self.center.same_level(point)
		xself = new_self[0]
		xpoint = new_point[0]
		cureps = eps / 2 * new_self[2]

		return xpoint == xself + cureps

	def get_top_neighbour(self):
		"""
		возвращает пиксель, смежный по верхней границе

		:return: экземпляр Pixel
		"""

		res = average(eps, self.center)
		new_eps = res[0]
		new_self = res[1]
		return get_pixel(Point(int(new_self.x), int(new_self.y + new_eps), int(new_self.z), homogeneous = True))

	def get_bottom_neighbour(self):
		"""
		возвращает пиксель, смежный по нижней границе

		:return: экземпляр Pixel
		"""

		res = average(eps, self.center)
		new_eps = res[0]
		new_self = res[1]
		return get_pixel(Point(int(new_self.x), int(new_self.y - new_eps), int(new_self.z), homogeneous = True))

	def get_left_neighbour(self):
		"""
		возвращает пиксель, смежный по левой границе

		:return: экземпляр Pixel
		"""

		res = average(eps, self.center)
		new_eps = res[0]
		new_self = res[1]
		return get_pixel(Point(int(new_self.x - new_eps), int(new_self.y), int(new_self.z), homogeneous = True))

	def get_right_neighbour(self):
		"""
		возвращает пиксель, смежный по правой границе

		:return: экземпляр Pixel
		"""

		res = average(eps, self.center)
		new_eps = res[0]
		new_self = res[1]
		return get_pixel(Point(int(new_self.x + new_eps), int(new_self.y), int(new_self.z), homogeneous = True))

	def get_neighbour(self, point):
		"""
		возвращает пиксель, смежный по той грани, которой принадлежит точка

		:param point: экземпляр Point
		:return: Pixel смежный пиксель
		"""	

		res = average(eps, self.center)
		new_eps = res[0]
		new_self = res[1]

		if self.is_on_top(point):
			return get_pixel(Point(int(new_self.x), int(new_self.y + new_eps), int(new_self.z), homogeneous = True))
		if self.is_on_bottom(point):
			return get_pixel(Point(int(new_self.x), int(new_self.y - new_eps), int(new_self.z), homogeneous = True))
		else:
			return None

	@property
	def x(self):
		return self.center.x

	@property
	def y(self):
		return self.center.y

	@property
	def z(self):
		return self.center.z

	@property
	def nw(self):
		res = average(eps / 2, self.center)
		new_eps = res[0]
		new_self = res[1]
		return Point(int(new_self.x - new_eps), int(new_self.y + new_eps), int(new_self.z), homogeneous = True)

	@property
	def ne(self):
		res = average(eps / 2, self.center)
		new_eps = res[0]
		new_self = res[1]
		return Point(int(new_self.x + new_eps), int(new_self.y + new_eps), int(new_self.z), homogeneous = True)

	@property
	def sw(self):
		res = average(eps / 2, self.center)
		new_eps = res[0]
		new_self = res[1]
		return Point(int(new_self.x - new_eps), int(new_self.y - new_eps), int(new_self.z), homogeneous = True)

	@property
	def se(self):
		res = average(eps / 2, self.center)
		new_eps = res[0]
		new_self = res[1]
		return Point(int(new_self.x + new_eps), int(new_self.y - new_eps), int(new_self.z), homogeneous = True)


class SweepLine:
	"""
	Класс заметающей прямой
	"""

	class Event:
		"""
		Класс события, обрабатываемого заметающей прямой
		"""

		class Type(Enum):
			"""
			Тип события
			"""

			SEG_END      = 0 # конец отрезка
			SEG_START    = 6 # начало отрезка
			SEG_SEG      = 2 # пересечение двух отрезков
			SEG_PIX      = 3 # пересечение отрезка с границей пикселя
			PIX_END      = 4 # конец пикселя
			SEG_REINSERT = 5 # повторное вхождение отрезка после выхода из пикселя

		def __init__(self, etype, point, segment=None, pixel=None, isstatus=None):
			"""
			Конструктор события

			:param etype: SweepLine.Event.Type тип события
			:param point: Point точка, в которой происходит событие
			:param segment: Segment для событий типа SEG_START, SEG_END, SEG_PIX, SEG_REINSERT \
			                -- отрезок, участвующий в них
			:param pixel: Pixel для событий типа SEG_PIX, PIX_END -- пиксель, участвующий в них
			:param isstatus: Происходит ли событие в массиве status или segments
			"""

			self.etype = etype
			self.point = normalize(point)
			self.segment = segment
			self.pixel = pixel
			self.isstatus = isstatus

		def __eq__(self, other):
			"""
			проверяет два события на равенство

			:param other: экземпляр SweepLine.Event 
			:return: bool
			"""

			return self.point == other.point and self.etype == other.etype \
			and self.segment == other.segment and self.pixel == other.pixel

		def __lt__(self, other):
			"""
			лексикографическое "меньше" для событий

			:param other: экземпляр SweepLine.Event
			:return: bool
			"""

			new_self, new_other = self.point.same_level(other.point)
			if self.point.x/self.point.z == other.point.x / other.point.z:
				if self.etype.value == other.etype.value:
					return self.point.y/self.point.z < other.point.y / other.point.z
				else:
					return self.etype.value < other.etype.value	
			return self.point.x/self.point.z < other.point.x / other.point.z

		def __str__(self):
			s = 'type: {}\npoint: {}'.format(self.etype.name, self.point)
			if self.segment is not None:
				s += '\nsegment: {}'.format(self.segment)
			if self.pixel is not None:
				s += '\npixel: {}'.format(self.pixel)
			return s

		def handle(self):
			"""
			вызывает функцию-обработчик для события данного типа
			"""

			handler = handlers[self.etype]
			if self.segment is not None and self.pixel is not None:
				handler(self.point, segment=self.segment, pixel=self.pixel)
			elif self.segment is not None:
				handler(self.point, segment=self.segment)
			elif self.pixel is not None:
				handler(self.point, pixel=self.pixel)
			else:
				handler(self.point, isstatus=self.isstatus)

		@property
		def x(self):
			return self.point.x

		@property
		def y(self):
			return self.point.y
		@property
		def z(self):
			return self.point.z


	def __init__(self, segments):
		"""
		Конструктор заметающей прямой

		:param segments: список отрезков Segment, для работы с которыми строится заметающая прямая
		"""

		#текущий x
		self.xpos = 0

		#массив вспомогательных отрезков, отсортированный по y координате при текущем x
		self.status = blist()

		#массив отрезков, пересекающих текущий x, отсортированный по y координате точки пересечения
		self.segments = blist()

		#массив событий
		self.events = []

		#массив пересечений отрезков из status, отсортированный по x
		self.intersections_status = SortedListWithKey(key = lambda xnpair: xnpair[0])
		
		#массив пересечений отрезков из segments, отсортированный по x
		self.intersections_segments = SortedListWithKey(key = lambda xnpair: xnpair[0])

		# инициализируем список событий началами и концами отрезков
		for s in segments:
			self.push(SweepLine.Event(SweepLine.Event.Type.SEG_START, s.start, segment=s))
			self.push(SweepLine.Event(SweepLine.Event.Type.SEG_END, s.end, segment=s))

	def bsearch(self, listt, segment):
		"""
		бинарный поиск, находит первый элемент больший данного

		:param listt: массив для поиска
		:param segment: экземпляр Segment
		:return: int
		"""

		r = len(listt) - 1
		if (r == -1):
			return 0
		searchX = segment.atX(self.xpos) + 0.00002
		hiddenX = searchX - 0.00004
		l = -1
		while (r - l > 1):
			m = (r + l) // 2
			if (listt[m].atX(self.xpos) > searchX):
				r = m
			else:
				l = m
		while (r in range(0, len(listt))) and (listt[r].atX(self.xpos) <= searchX):
			r += 1
		q = r 
		while ((q - 1) in range(0, len(listt))) and (listt[q - 1].atX(self.xpos) >= hiddenX) and ((listt[q - 1] != segment) or (listt[q - 1].isbound != segment.isbound)):
			q -= 1
		if ((q - 1) in range(0, len(listt))) and (listt[q - 1] == segment and listt[q - 1].isbound == segment.isbound):
			return q
		return r

	def sort_intersection(self, listt, deque_inter):
		"""
		обрабатывает пересечения отрезков для поддержания сортированности

		:param listt: массив отрезков
		:param deque_inter: массив пересечений
		"""

		pair = deque_inter.pop(0)
		a = pair[1]
		b = pair[2]
		ia = self.bsearch(listt, a) - 1
		ib = self.bsearch(listt, b) - 1
		position = min(a.end.x / a.end.z, b.end.x / b.end.z)
		if (ia < ib and a.atX(position) > b.atX(position)) or (ia > ib and a.atX(position) < b.atX(position)):
			listt.pop(ia)
			listt.insert(ia, b)
			listt.pop(ib)
			listt.insert(ib, a)

	def find_intersections(self, listt, segment, start, delta, deque_inter, isstatus, msg):
		"""
		находит все пересечения отрезка segment с отрезками в массиве listt начиная с индекса start и шагом delta

		:param listt: массив отрезков
		:param segment: экземпляр Segment
		:param start: int
		:param delta: int
		:param deque_inter: массив пересечений
		:param isstatus: bool
		:param msg: string
		"""

		j = start
		while (j + delta in range(0, len(listt))):
			p = segment.intersects(listt[j + delta])
			if p is not None and (p.x / p.z >= self.xpos):
				self.push(SweepLine.Event(SweepLine.Event.Type.SEG_SEG, p, isstatus=isstatus))
				deque_inter.add((p, segment, listt[j + delta]))
				j += delta
			else:
				break

	def insert(self, listt, segment, deque_inter, isstatus, msg):
		"""
		добавляет отрезок segment в массив listt, также вычисляет новые пересечения и добавляет их в очередь deque_inter

		:param listt: массив
		:param segment: экземпляр Segment
		:param deque_inter: массив 
		:param isstatus: bool
		:param msg: string
		"""

		i = self.bsearch(listt, segment)
		if (i in range(1, len(listt) + 1) and listt[i - 1] == segment):
			i -= 1
		else:
			listt.insert(i, segment)
		self.find_intersections(listt, segment, i, -1, deque_inter, isstatus, msg)
		self.find_intersections(listt, segment, i, +1, deque_inter, isstatus, msg)

	def remove(self, listt, segment, deque_inter, isstatus, msg):
		"""
		удаляет отрезок segment из массива listt, также вычисляет новые пересечения и добавляет их в очередь deque_inter

		:param listt: массив
		:param segment: экземпляр Segment
		:param deque_inter: массив 
		:param isstatus: bool
		:param msg: string
		"""

		i = self.bsearch(listt, segment) - 1
		while (i in range (0, len(listt)) and listt[i] == segment and listt[i].isbound == segment.isbound):
			listt.pop(i)
			i -= 1
		if (i == -1):
			return
		self.find_intersections(listt, listt[i], i, +1, deque_inter, isstatus, msg)
		if (i == len(listt) - 1):
			return
		self.find_intersections(listt, listt[i + 1], i + 1, -1, deque_inter, isstatus, msg)


	def push(self, event):
		"""
		добавляет событие в очередь events

		:param event: экземпляр SweepLine.Event
		"""

		heapq.heappush(self.events, event)

	def peek(self):
		"""
		возвращает самое раннее событие из events, не удаляя его из очереди

		:return: SweepLine.Event
		"""

		return self.events[0]

	def pop(self):
		"""
		возвращает самое раннее событие из events, удаляя его из очереди

		:return: SweepLine.Event
		"""

		return heapq.heappop(self.events)

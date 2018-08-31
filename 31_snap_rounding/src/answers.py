from sortedcontainers import SortedListWithKey
from structs import Segment, Pixel, SweepLine as SL, Point, normalize, get_pixel, add_pixel_to_seg
from testing import hot, current, line, starting_segments
from collections import deque


# возвращает все прилежащие горячие пиксели если точка на верхней/нижней границе
def bcheck(point):
	pixel = Pixel(point)
	ans = []
	if (pixel.is_on_top(point) and pixel in current):
		ans.append(pixel)
	if (pixel.is_on_top(point) and pixel.get_top_neighbour() in current):
		ans.append(pixel.get_top_neighbour())
	return ans



def segment_endpoint_answer(point, segment):
	pixel = get_pixel(point)
	if pixel not in current:
		# до этого в пикселе не было обнаружено критических точек
		heat_answer(pixel)

	# если событие -- начало отрезка
	if point == segment.start:
		line.insert(line.segments, segment, line.intersections_segments, isstatus=False, msg="segment")
		add_pixel_to_seg(pixel, segment)
		pixel.segs.append(segment)

		# игнорируем ту часть отрезка, что лежит внутри пикселя
		for intersection in segment.intersections(pixel):
			if intersection is not None:
				line.push(SL.Event(SL.Event.Type.SEG_REINSERT, intersection, segment=segment, pixel=pixel))
				break;

	# если событие -- конец отрезка
	else:
		line.remove(line.status, segment, line.intersections_status, isstatus=True, msg="status")
		line.remove(line.segments, segment, line.intersections_segments, isstatus=False, msg="segment")

def segseg_intersection_answer(point, isstatus):

	# если пересечение отрезков из segments
	if (not isstatus): 
		line.sort_intersection(line.segments, line.intersections_segments)
		return
	pair = line.intersections_status[0]
	line.sort_intersection(line.status, line.intersections_status)

	# рассматриваем случай пересечения обычного отрезка и границы пикселя
	if (pair[1].isbound != 0 or pair[2].isbound != 0):
		return
	new_point = normalize(point)
	pixel = get_pixel(new_point)
	if pixel not in current:
		heat_answer(pixel)

def segpix_intersection_answer(point, segment, pixel):
	if pixel.is_on_top(point):
		pixel.upper.append(segment)
	else:
		pixel.lower.append(segment)

	# добавляем отрезок в соответствующий список
	add_pixel_to_seg(pixel, segment)
	pixel.segs.append(segment)

	# игнорируем фрагмент отрезка, находящийся внутри пикселя
	line.remove(line.status, segment, line.intersections_status, isstatus=True, msg="status")
	
	# ищем пересечение с пикселем
	lastinter = None
	for intersection in segment.intersections(pixel):
		if ((intersection is not None) and (intersection > point)) and ((lastinter is None) or (intersection > lastinter)):
			lastinter = intersection
	if (lastinter is not None):
		line.push(SL.Event(SL.Event.Type.SEG_REINSERT, lastinter, segment=segment, pixel=pixel))


def segment_reinsertion_answer(point, segment, pixel):
	line.insert(line.segments, segment, line.intersections_segments, isstatus=False, msg="segment")
	line.insert(line.status, segment, line.intersections_status, isstatus=True, msg="status")
	if pixel.is_on_top(point):
		pixel.upper.append(segment)
	else:
		pixel.lower.append(segment)
	neighbour = pixel.get_neighbour(point)
	if (neighbour is not None) and (neighbour in current):
		segpix_intersection_answer(point, segment, neighbour)

def pixel_end_answer(point, pixel):

	# удаляем границы пикселя из status
	if pixel.is_on_top(point):
		line.remove(line.status, pixel.top(), line.intersections_status, isstatus=True, msg="status")
	else:
		line.remove(line.status, pixel.bottom(), line.intersections_status, isstatus=True, msg="status")
	hot.extend(current)
	current.clear()

def heat_answer(pixel):
	pixel.center = normalize(pixel.center)
	l = []
	maybegood = []

	# найдем горячий пиксель выше и ниже текущего
	mypos = current.bisect(pixel)
	lower = None
	higher = None
	if not (0 == len(current)):
		if (mypos > 0 and mypos < len(current)):
			lower = current[mypos - 1]
			maybegood.extend(lower.upper)
			higher = current[mypos]
			maybegood.extend(higher.lower)
		if (mypos == 0):
			higher = current[mypos]
			maybegood.extend(higher.lower)
		else:
			if (mypos == len(current)):
				lower = current[len(current) - 1]
				maybegood.extend(lower.upper)
	
	current.add(pixel)
	
	# найдем такой отрезок [li; hi] (индексы массива segments, отсортированного по y) что значения элементов находятся
	# между нижней границей верхнего пикселя и верхней границей нижнего
	li = 0
	if (lower is not None):
		li = max(li, line.bsearch(line.segments, lower.top()))
		while (li in range(1, len(line.segments)) and line.segments[li].atX(line.xpos) >= lower.top().atX(line.xpos)):
			li -= 1
	hi = len(line.segments)
	if (higher is not None):
		hi = min(hi, line.bsearch(line.segments, higher.bottom()))
		while (hi in range(0, len(line.segments)) and line.segments[hi].atX(line.xpos) <= higher.bottom().atX(line.xpos)):
			hi += 1
	extender = line.segments[li:hi]
	maybegood.extend(extender)

	# из всех возможных отрезков берем те, что пересекаются с пикселем и запоминаем это
	for s in maybegood:
		if (s.isbound != 0):
			continue;
		for intersection in s.intersections(pixel) :
			if (intersection is not None):
				add_pixel_to_seg(pixel, s)
				pixel.segs.append(s)
				l.append(s)
				break

	# из всех таких берем те, что имеют пересечение после xpos и удаляем фрагмент внутри пикселя
	for segment in l:
		somearr = []
		for intersection in segment.intersections(pixel):
			if intersection is not None and (intersection.x / intersection.z) > line.xpos:
				somearr.append(intersection)
		if len(somearr) >= 2:
			if (somearr[0] < somearr[1]):
				swap(somearr[0], somearr[1])
		if len(somearr) >= 1:
			line.remove(line.status, segment, line.intersections_status, isstatus=True, msg="status")
			line.push(SL.Event(SL.Event.Type.SEG_REINSERT, somearr[0], segment=segment, pixel=pixel))

	# собираем bottom и top нового пикселя
	for segment in pixel.segs:
		inter = segment.intersects(pixel.bottom())
		if inter is not None:
			pixel.lower.append(segment)
		inter = segment.intersects(pixel.top())
		if inter is not None:
			pixel.upper.append(segment)

	# добавляем границы пикселя в status
	top, bottom = pixel.top(), pixel.bottom()
	line.insert(line.status, top, line.intersections_status, isstatus=True, msg="status")
	line.push(SL.Event(SL.Event.Type.PIX_END, top.end, pixel=pixel))
	line.insert(line.status, bottom, line.intersections_status, isstatus=True, msg="status")
	line.push(SL.Event(SL.Event.Type.PIX_END, bottom.end, pixel=pixel))

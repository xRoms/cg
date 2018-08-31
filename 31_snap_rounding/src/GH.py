from test_parameters import handlers
from structs import SweepLine as SL
from testing import line
from answers import *

def segment_endpoint(point, segment):
	segment_endpoint_answer(point, segment)

def segseg_intersection(point, isstatus):
	segseg_intersection_answer(point, isstatus)

def segpix_intersection(point, segment, pixel):
	segpix_intersection_answer(point, segment, pixel)

def segment_reinsertion(point, segment, pixel):
	segment_reinsertion_answer(point, segment, pixel)

def pixel_end(point, pixel):
	pixel_end_answer(point, pixel)

handlers[SL.Event.Type.SEG_END] = segment_endpoint
handlers[SL.Event.Type.SEG_START] = segment_endpoint
handlers[SL.Event.Type.SEG_SEG] = segseg_intersection
handlers[SL.Event.Type.SEG_PIX] = segpix_intersection
handlers[SL.Event.Type.SEG_REINSERT] = segment_reinsertion
handlers[SL.Event.Type.PIX_END] = pixel_end

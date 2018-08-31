import re
import os
import shutil
import struct
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import collections as mc
import ipywidgets as widgets
from collections import defaultdict

from test_parameters import eps, bound
from structs import SweepLine as SL, Segment, Point, pixelspassed
from testing import line, hot, current

def _type_color(e):
	if e.etype == SL.Event.Type.SEG_END:
		return '>g'
	if e.etype == SL.Event.Type.SEG_START:
		return '<g'
	if e.etype == SL.Event.Type.SEG_SEG:
		return 'or'
	if e.etype == SL.Event.Type.SEG_PIX:
		return 'хb'
	if e.etype == SL.Event.Type.PIX_END:
		return 'xb'
	if e.etype == SL.Event.Type.SEG_REINSERT:
		return '<r'

def draw_all(segments, status, events, hot, current, yasegments):
	ax_min, ax_max = 0 - eps / 2, bound + eps / 2 
	fig = plt.figure(figsize=(8,8))
	ax = fig.add_subplot(111, aspect='equal')
	ax.set_xlim(ax_min, ax_max)
	ax.set_ylim(ax_min, ax_max)
	ax.set_xticks(np.arange(-eps / 2, bound, eps), minor=True)
	ax.set_yticks(np.arange(-eps / 2, bound, eps), minor=True)
	ax.grid(which='minor')
	ax.grid(which='major', alpha=0.0)

	seg_col = []
	for s in status:
		seg_col.append('green')
	seg_src = list(map(lambda seg: [(seg.start.x / seg.start.z, seg.start.y / seg.start.z), (seg.end.x / seg.end.z, seg.end.y / seg.end.z)], status))
	ax.add_collection(mc.LineCollection(seg_src, colors=seg_col))

	seg_col = []
	seg_src = []
	for s in segments:
		if s not in status:
			seg_col.append('blue')
			seg_src.append(s)
	seg_src = list(map(lambda seg: [(seg.start.x / seg.start.z, seg.start.y / seg.start.z), (seg.end.x / seg.end.z, seg.end.y / seg.end.z)], seg_src))
	ax.add_collection(mc.LineCollection(seg_src, colors=seg_col))


	for p in current:
		ax.add_patch(patches.Rectangle((p.sw.x / p.sw.z, p.sw.y / p.sw.z), eps, eps, alpha=0.5, facecolor="yellow"))

	for p in hot:
		ax.add_patch(patches.Rectangle((p.sw.x / p.sw.z, p.sw.y / p.sw.z), eps, eps, alpha=0.5, facecolor="red"))

	e_src = defaultdict(list)
	for e in events:
		e_src[_type_color(e)].append(e)
	for t in e_src:
		plt.plot([(e.x / e.z) for e in e_src[t]], [(e.y / e.z) for e in e_src[t]], t, label=e_src[t][0].etype.name, zorder=1)
	plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,numpoints=1)


	# plt.show()
	return fig

def visual_dump_pieces(xpos, status, segments, hot, current, yasegments, filename):	
	fig = draw_all(segments, status, line.events, hot, current, yasegments)
	if xpos is not None:
		plt.axvline(xpos, ymin=0.05, ymax=0.95, color='black', linestyle='dashed')
	plt.savefig(filename)
	plt.close(fig)

def dump_answer(segments, hot, filename):	
	fig = draw_result(segments, hot)
	plt.savefig(filename)
	plt.close(fig)

def create_dump_func(folder, func, *args):
	c = 0
	shutil.rmtree(folder, ignore_errors=True)
	os.makedirs(folder, exist_ok=True)
	def dump(x=None):
		nonlocal c
		func(x, *args, filename = os.path.join(folder,'{}.png'.format(c)))
		c += 1
	return dump

def create_dump_answer(folder, func, *args):
	c = 0
	shutil.rmtree(folder, ignore_errors=True)
	os.makedirs(folder, exist_ok=True)
	def dump():
		nonlocal c
		func(*args, filename = os.path.join(folder,'{}.png'.format(c)))
		c += 1
	return dump

def natural_sort(l): 
	convert = lambda text: int(text) if text.isdigit() else text.lower() 
	alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
	return sorted(l, key = alphanum_key)

def _get_png_info(data):
	w, h = struct.unpack('>LL', data[16:24])
	width = int(w)
	height = int(h)
	return width, height
    
def SlideShower(folder, frame_duration=800):
	slides = list([open(os.path.join(folder,s), 'rb').read()
                   for s in natural_sort(os.listdir(folder))])

	x, y = _get_png_info(slides[0])
	img = widgets.Image(value=slides[0], width=x, height=y)

	def on_frame(change):
		n = change['new']
		img.value = slides[n]

	play = widgets.Play(
        value=0,
        min=0,
        max=len(slides),
        step=1,
        interval=frame_duration,
        disabled=False
    )
	slider = widgets.IntSlider(
        value=0,
        min=0,
        max=len(slides)-1,
        step=1
    )
	slider.observe(on_frame, names='value')
	widgets.jslink((play, 'value'), (slider, 'value'))
	box = widgets.VBox([img, widgets.HBox([play, slider])])
	return box

def draw_result(segments, hot):
	ax_min, ax_max = 0 - eps / 2, bound + eps / 2 
	fig = plt.figure(figsize=(8,8))
	ax = fig.add_subplot(111, aspect='equal')
	ax.set_xlim(ax_min, ax_max)
	ax.set_ylim(ax_min, ax_max)
	ax.set_xticks(np.arange(-eps / 2, bound, eps), minor=True)
	ax.set_yticks(np.arange(-eps / 2, bound, eps), minor=True)
	ax.grid(which='minor')
	ax.grid(which='major', alpha=0.0)

	for p in hot:
		ax.add_patch(patches.Rectangle((p.sw.x / p.sw.z, p.sw.y / p.sw.z), eps, eps, alpha=0.5, facecolor="red"))

	seg_col = []
	seg_src = []
	for p in pixelspassed.keys():
		lis = pixelspassed[p]
		lis.sort(key=lambda smth: smth[0])
		q = 0
		while (q < len(lis) - 1):
			seg = Segment(Point(int(lis[q][1][0]), int(lis[q][1][1]), int(lis[q][1][2]), homogeneous=True), Point(int(lis[q + 1][1][0]), int(lis[q + 1][1][1]), int(lis[q + 1][1][2]), homogeneous=True))
			seg_col.append('black')
			seg_src.append(seg)
			q += 1

	seg_src = list(map(lambda seg: [(seg.start.x / seg.start.z, seg.start.y / seg.start.z), (seg.end.x / seg.end.z, seg.end.y / seg.end.z)], seg_src))
	ax.add_collection(mc.LineCollection(seg_src, colors=seg_col))		

	return fig


def Result(folder):
	this = open(os.path.join(folder,"0.png"), 'rb').read()
	x, y = _get_png_info(this)
	img = widgets.Image(value=this, width=x, height=y)
	return img
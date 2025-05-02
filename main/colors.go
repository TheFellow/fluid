package main

import (
	"image/color"
	"math"
)

func getSciValue(val, minVal, maxVal float32) color.RGBA {
	val = min(max(val, minVal), maxVal-0.0001)
	var d = maxVal - minVal
	if d <= 0 {
		val = 0.5
	} else {
		val = (val - minVal) / d
	}
	var m = float32(0.25)
	var num = float32(math.Floor(float64(val / m)))
	var s = (val - num*m) / m
	var r, g, b float32

	switch num {
	case 0:
		r = 0.0
		g = s
		b = 1.0
	case 1:
		r = 0.0
		g = 1.0
		b = 1.0 - s
	case 2:
		r = s
		g = 1.0
		b = 0.0
	case 3:
		r = 1.0
		g = 1.0 - s
		b = 0.0
	}

	return color.RGBA{
		R: uint8(255 * r),
		G: uint8(255 * g),
		B: uint8(255 * b),
		A: 0xff,
	}
}

package main

import (
	"image/color"
	"math"
)

// getDivergingColor maps a value to a blue-white-red diverging colormap.
// Negative values are blue, zero is white, positive values are red.
func getDivergingColor(val, minVal, maxVal float32) color.RGBA {
	// Normalize to [-1, 1] where 0 is the midpoint
	absMax := float32(math.Max(math.Abs(float64(minVal)), math.Abs(float64(maxVal))))
	if absMax < 1e-8 {
		return color.RGBA{R: 255, G: 255, B: 255, A: 255}
	}
	t := val / absMax // [-1, 1]
	if t > 1 {
		t = 1
	}
	if t < -1 {
		t = -1
	}

	var r, g, b float32
	if t >= 0 {
		// White to red
		r = 1.0
		g = 1.0 - t
		b = 1.0 - t
	} else {
		// White to blue
		a := -t
		r = 1.0 - a
		g = 1.0 - a
		b = 1.0
	}

	return color.RGBA{
		R: uint8(255 * r),
		G: uint8(255 * g),
		B: uint8(255 * b),
		A: 0xff,
	}
}

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

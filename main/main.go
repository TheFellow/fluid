package main

import (
	"fmt"
	"image"
	"log"
	"math"

	"github.com/TheFellow/fluid/pkg/fluid"
	"github.com/hajimehoshi/ebiten/v2"
	"github.com/hajimehoshi/ebiten/v2/ebitenutil"
)

const (
	screenWidth  = 1024
	screenHeight = 768
	fluidWidth   = 400
	fluidHeight  = 280
)

type Game struct {
	fluid *fluid.Fluid
	image *image.RGBA
}

func NewGame() *Game {
	f := fluid.New(1000, fluidWidth, fluidHeight, 1.0/100.0)
	// Set no walls except floor
	for i := range fluidWidth {
		for j := range fluidHeight {
			if i == 0 || i == fluidWidth-1 || j == 0 || j == fluidHeight-1 {
				f.SetSolid(i, j, true)
			} else {
				f.SetSolid(i, j, false)
			}
		}
	}
	return &Game{
		fluid: f,
		image: image.NewRGBA(image.Rect(0, 0, fluidWidth+2, fluidHeight+2)),
	}
}

func (g *Game) Update() error {
	dt := float32(1.0 / float32(120))
	g.fluid.Simulate(dt, -9.81, 10)
	return nil
}

var drawOpts = &ebiten.DrawImageOptions{Filter: ebiten.FilterLinear}

func (g *Game) Draw(screen *ebiten.Image) {
	// min/max pressures
	minPressure := float32(math.MaxFloat32)
	maxPressure := float32(0)
	pressures := g.fluid.Pressure()
	for i := range fluidWidth {
		for j := range fluidHeight {
			p, err := pressures.Value(i, j)
			if err != nil {
				log.Panicf("cannot get pressure: %v", err)
			}
			if p < minPressure {
				minPressure = p
			}
			if p > maxPressure {
				maxPressure = p
			}
		}
	}

	for i := range pressures.NumX {
		for j := range pressures.NumY {
			p, err := pressures.Value(i, fluidHeight+2-j-1)
			if err != nil {
				log.Panicf("cannot get pressure: %v", err)
			}
			color := getSciValue(p, minPressure, maxPressure)
			g.image.Pix[4*i+j*g.image.Stride+0] = color.R
			g.image.Pix[4*i+j*g.image.Stride+1] = color.G
			g.image.Pix[4*i+j*g.image.Stride+2] = color.B
			g.image.Pix[4*i+j*g.image.Stride+3] = color.A
		}
	}

	// v := g.fluid.Velocity()
	// x, y, err := v.Value(fluidWidth/2, fluidHeight/2)
	// if err != nil {
	// 	log.Panicf("cannot get velocity: %v", err)
	// }

	// render
	eImage := ebiten.NewImageFromImage(g.image)
	//screen.WritePixels(g.image.Pix)
	screen.DrawImage(eImage, drawOpts)
	ebitenutil.DebugPrint(screen,
		fmt.Sprintf("FluidSim - FPS: %0.2f", ebiten.ActualFPS()),
	)
}

func (g *Game) Layout(outsideWidth, outsideHeight int) (w, hHeight int) {
	return g.fluid.NumX, g.fluid.NumY
}

func main() {
	ebiten.SetWindowSize(screenWidth, screenHeight)
	ebiten.SetWindowTitle("FluidSim")

	if err := ebiten.RunGame(NewGame()); err != nil {
		log.Fatal(err)
	}
}

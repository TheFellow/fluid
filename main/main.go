package main

import (
	"fmt"
	"log"
	"math"

	"github.com/TheFellow/fluid/pkg/fluid"
	"github.com/hajimehoshi/ebiten/v2"
	"github.com/hajimehoshi/ebiten/v2/ebitenutil"
)

const (
	screenWidth  = 800
	screenHeight = 480
	fluidWidth   = 400
	fluidHeight  = 240
)

type Game struct {
	fluid *fluid.Fluid
}

func NewGame() *Game {
	f := fluid.New(1000, fluidWidth, fluidHeight, float32(screenHeight)/200.0)
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
	}
}

func (g *Game) Update() error {
	dt := float32(1.0 / float32(ebiten.TPS()))
	g.fluid.Simulate(dt, -9.81, 10)
	return nil
}

func (g *Game) Draw(screen *ebiten.Image) {
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
	for i := range fluidWidth {
		for j := range fluidHeight {
			_, _ = i, j
		}
	}

	v := g.fluid.Velocity()
	x, y, err := v.Value(fluidWidth/2, fluidHeight/2)
	if err != nil {
		log.Panicf("cannot get velocity: %v", err)
	}

	ebitenutil.DebugPrint(screen,
		fmt.Sprintf("FluidSim - FPS: %0.2f\nPressures: [%0.2f, %0.2f]\nVelocity: (%0.2f, %0.2f)",
			ebiten.ActualFPS(), minPressure, maxPressure, x, y),
	)
}

func (g *Game) Layout(outsideWidth, outsideHeight int) (w, hHeight int) {
	return fluidWidth, fluidHeight
}

func main() {
	ebiten.SetWindowSize(screenWidth, screenHeight)
	ebiten.SetWindowTitle("FluidSim")

	if err := ebiten.RunGame(NewGame()); err != nil {
		log.Fatal(err)
	}
}

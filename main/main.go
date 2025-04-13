package main

import (
	"fmt"
	"image/color"
	"log"
	"math/rand/v2"

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
	counter int
	bgColor color.RGBA

	fluid *fluid.Fluid
}

func NewGame() *Game {
	f := fluid.New(1000, fluidWidth, fluidHeight, 1000)
	return &Game{
		fluid: f,
	}
}

func (g *Game) Update() error {
	g.counter++
	if g.counter%ebiten.TPS() == 0 {
		g.bgColor = color.RGBA{
			R: uint8(rand.IntN(256)),
			G: uint8(rand.IntN(256)),
			B: uint8(rand.IntN(256)),
			A: 0xff,
		}
	}

	dt := float32(1.0 / ebiten.ActualTPS())
	g.fluid.Simulate(dt, -9.81, 10)
	return nil
}

func (g *Game) Draw(screen *ebiten.Image) {
	screen.Fill(g.bgColor)

	ebitenutil.DebugPrint(screen, fmt.Sprintf("FluidSim - FPS: %0.2f", ebiten.ActualFPS()))
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

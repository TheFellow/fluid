package main

import (
	"fmt"
	"image/color"
	"log"
	"math/rand/v2"
	"time"

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
}

func NewGame() *Game {
	return &Game{}
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
	return nil
}

func (g *Game) Draw(screen *ebiten.Image) {
	screen.Fill(g.bgColor)

	ebitenutil.DebugPrint(screen, fmt.Sprintf("FluidSim - %s\nFPS: %0.2f",
		time.Now().Format("2006-01-02 15:04:05"), ebiten.ActualFPS()))
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

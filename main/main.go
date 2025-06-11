package main

import (
	"fmt"
	"image"
	"log"

	"github.com/TheFellow/fluid/pkg/fluid"
	"github.com/hajimehoshi/ebiten/v2"
	"github.com/hajimehoshi/ebiten/v2/ebitenutil"
	"github.com/hajimehoshi/ebiten/v2/inpututil"
)

const (
	factor = 4

	fluidWidth   = 300
	fluidHeight  = 251
	screenWidth  = factor * fluidWidth
	screenHeight = factor * fluidHeight
)

type Game struct {
	fluid *fluid.Fluid
	image *image.RGBA

	showSmoke bool
	jet       bool

	dragging  bool
	dragValue bool
	prevI     int
	prevJ     int
}

func NewGame() *Game {
	f := fluid.New(1000, fluidWidth, fluidHeight, 1.0/100.0)
	// Configure walls
	for i := range fluidWidth + 2 {
		for j := range fluidHeight + 2 {
			if j == 0 || j == fluidHeight+2-1 || i == 0 /* || i == fluidWidth+2-1 */ {
				f.SetSolid(i, j, true)
			} else {
				f.SetSolid(i, j, false)
			}
		}
	}
	return &Game{
		fluid:     f,
		image:     image.NewRGBA(image.Rect(0, 0, fluidWidth+2, fluidHeight+2)),
		showSmoke: true,
	}
}

func (g *Game) Update() error {
	if inpututil.IsKeyJustPressed(ebiten.KeyS) {
		g.showSmoke = !g.showSmoke
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyJ) {
		g.jet = !g.jet
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyG) {
		if g.fluid.Gravity != 0.0 {
			g.fluid.SetGravity(false)
		} else {
			g.fluid.SetGravity(true)
		}
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyR) {
		g.fluid.Reset()
	}

	if inpututil.IsMouseButtonJustPressed(ebiten.MouseButtonLeft) {
		mx, my := ebiten.CursorPosition()
		i, j := g.screenToFluid(mx, my)
		if g.inBounds(i, j) {
			g.dragging = true
			g.dragValue = !g.fluid.IsSolid(i, j)
			g.drawLine(i, j, i, j, g.dragValue)
			g.prevI = i
			g.prevJ = j
		}
	} else if ebiten.IsMouseButtonPressed(ebiten.MouseButtonLeft) && g.dragging {
		mx, my := ebiten.CursorPosition()
		i, j := g.screenToFluid(mx, my)
		if g.inBounds(i, j) && (i != g.prevI || j != g.prevJ) {
			g.drawLine(g.prevI, g.prevJ, i, j, g.dragValue)
			g.prevI = i
			g.prevJ = j
		}
	} else if g.dragging {
		g.dragging = false
	}

	if g.jet {
		x, y := 4.0, 0.0
		size := 100
		for j := fluidHeight/2 - size; j < fluidHeight/2+size; j++ {
			g.fluid.SetVelocity(1, j, float32(x), float32(y))
			g.fluid.AddSmoke(1, j, 1.0)
		}
	}

	dt := float32(1.0 / float32(120))
	g.fluid.Simulate(dt, g.fluid.Gravity, 20)
	return nil
}

var drawOpts = &ebiten.DrawImageOptions{Filter: ebiten.FilterLinear}

func (g *Game) Draw(screen *ebiten.Image) {

	if !g.showSmoke {
		pressures := g.fluid.Pressure()
		for i := range pressures.NumX {
			for j := range pressures.NumY {
				p, err := pressures.Value(i, pressures.NumY-j-1)
				if err != nil {
					log.Panicf("cannot get pressure: %v", err)
				}
				color := getSciValue(p, pressures.MinValue, pressures.MaxValue)
				idx := g.fluidToImageIndex(i, j)
				g.image.Pix[idx+0] = color.R
				g.image.Pix[idx+1] = color.G
				g.image.Pix[idx+2] = color.B
				g.image.Pix[idx+3] = color.A
			}
		}
	}

	if g.showSmoke {
		smoke := g.fluid.Smoke()

		for i := range smoke.NumX {
			for j := range smoke.NumY {
				s, err := smoke.Value(i, smoke.NumY-j-1)
				if err != nil {
					log.Panicf("cannot get smoke: %v", err)
				}
				color := getSciValue(s, smoke.MinValue, smoke.MaxValue)
				idx := g.fluidToImageIndex(i, j)
				g.image.Pix[idx+0] = color.R
				g.image.Pix[idx+1] = color.G
				g.image.Pix[idx+2] = color.B
				g.image.Pix[idx+3] = color.A
			}
		}
	}

	for i := range g.fluid.NumX {
		for j := range g.fluid.NumY {
			if g.fluid.IsSolid(i, g.fluid.NumY-j-1) {
				idx := g.fluidToImageIndex(i, j)
				g.image.Pix[idx+0] = 0
				g.image.Pix[idx+1] = 0
				g.image.Pix[idx+2] = 0
				g.image.Pix[idx+3] = 0xff
			}
		}
	}

	// render
	eImage := ebiten.NewImageFromImage(g.image)
	screen.DrawImage(eImage, drawOpts)
	show := "Pressure"
	if g.showSmoke {
		show = "Smoke"
	}
	ebitenutil.DebugPrint(screen,
		fmt.Sprintf("FPS: %0.2f [S]how: %v\n[J]et: %v [G]ravity: %0.2f",
			ebiten.ActualFPS(), show, g.jet, g.fluid.Gravity),
	)
}

func (g *Game) fluidToImageIndex(i, j int) int {
	return 4*i + j*g.image.Stride
}

func (g *Game) screenToFluid(mx, my int) (int, int) {
	// CursorPosition already returns coordinates in the game's logical
	// resolution defined by Layout. As Layout returns the fluid grid size,
	// simply flip the y axis to convert to the fluid coordinates.
	i := mx
	j := g.fluid.NumY - my - 1
	return i, j
}

func (g *Game) inBounds(i, j int) bool {
	return i >= 0 && i < g.fluid.NumX && j >= 0 && j < g.fluid.NumY
}

func (g *Game) drawLine(x0, y0, x1, y1 int, value bool) {
	dx := abs(x1 - x0)
	dy := -abs(y1 - y0)
	sx := -1
	if x0 < x1 {
		sx = 1
	}
	sy := -1
	if y0 < y1 {
		sy = 1
	}
	err := dx + dy
	for {
		if g.inBounds(x0, y0) {
			g.fluid.SetSolid(x0, y0, value)
		}
		if x0 == x1 && y0 == y1 {
			break
		}
		e2 := 2 * err
		if e2 >= dy {
			err += dy
			x0 += sx
		}
		if e2 <= dx {
			err += dx
			y0 += sy
		}
	}
}

func abs(a int) int {
	if a < 0 {
		return -a
	}
	return a
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

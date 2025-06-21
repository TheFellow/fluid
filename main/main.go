package main

import (
	"fmt"
	"image"
	"image/color"
	"log"
	"runtime"
	"sync"

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

// parallelRange executes fn for each i in [start,end). The range is split among
// available CPUs.
func parallelRange(start, end int, fn func(i int)) {
	total := end - start
	if total <= 0 {
		return
	}
	workers := runtime.GOMAXPROCS(0)
	if workers > total {
		workers = total
	}
	var wg sync.WaitGroup
	chunk := (total + workers - 1) / workers
	for w := 0; w < workers; w++ {
		s := start + w*chunk
		e := s + chunk
		if e > end {
			e = end
		}
		if s >= end {
			break
		}
		wg.Add(1)
		go func(ss, ee int) {
			for i := ss; i < ee; i++ {
				fn(i)
			}
			wg.Done()
		}(s, e)
	}
	wg.Wait()
}

type Game struct {
	fluid  *fluid.Fluid
	image  *image.RGBA
	eImage *ebiten.Image

	showSmoke bool
	jet       bool

	showArrows bool
	paused     bool
	speed      float32

	wallTop    bool
	wallBottom bool
	wallLeft   bool
	wallRight  bool

	dragging  bool
	dragValue bool
	prevI     int
	prevJ     int
}

func NewGame() *Game {
	f := fluid.New(1000, fluidWidth, fluidHeight, 1.0/100.0)
	img := image.NewRGBA(image.Rect(0, 0, fluidWidth+2, fluidHeight+2))
	g := &Game{
		fluid:      f,
		image:      img,
		eImage:     ebiten.NewImageFromImage(img),
		showSmoke:  true,
		speed:      1.0,
		wallTop:    true,
		wallBottom: true,
		wallLeft:   true,
		wallRight:  false,
	}

	// Initialize all cells as empty and apply the boundary walls
	for i := 0; i < g.fluid.NumX; i++ {
		for j := 0; j < g.fluid.NumY; j++ {
			g.fluid.SetSolid(i, j, false)
		}
	}
	g.applyWallSettings()
	return g
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
	if inpututil.IsKeyJustPressed(ebiten.KeyV) {
		if g.fluid.Confinement != 0 {
			g.fluid.Confinement = 0
		} else {
			g.fluid.Confinement = 5
		}
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyA) {
		g.showArrows = !g.showArrows
	}
	if inpututil.IsKeyJustPressed(ebiten.KeySpace) {
		g.paused = !g.paused
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyEqual) || inpututil.IsKeyJustPressed(ebiten.KeyKPAdd) {
		g.speed *= 1.1
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyMinus) || inpututil.IsKeyJustPressed(ebiten.KeyKPSubtract) {
		g.speed /= 1.1
		if g.speed < 0.01 {
			g.speed = 0.01
		}
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyR) {
		g.fluid.Reset()
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyUp) {
		g.wallTop = !g.wallTop
		g.applyWallSettings()
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyDown) {
		g.wallBottom = !g.wallBottom
		g.applyWallSettings()
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyLeft) {
		g.wallLeft = !g.wallLeft
		g.applyWallSettings()
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyRight) {
		g.wallRight = !g.wallRight
		g.applyWallSettings()
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyC) {
		g.clearWalls()
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

	dt := float32(1.0/120.0) * g.speed
	if !g.paused {
		g.fluid.Simulate(dt, 20)
	}
	return nil
}

var drawOpts = &ebiten.DrawImageOptions{Filter: ebiten.FilterLinear}

func (g *Game) Draw(screen *ebiten.Image) {

	if !g.showSmoke {
		pressures := g.fluid.Pressure()
		parallelRange(0, pressures.NumX, func(i int) {
			for j := 0; j < pressures.NumY; j++ {
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
		})
	}

	if g.showSmoke {
		smoke := g.fluid.Smoke()

		parallelRange(0, smoke.NumX, func(i int) {
			for j := 0; j < smoke.NumY; j++ {
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
		})
	}

	parallelRange(0, g.fluid.NumX, func(i int) {
		for j := 0; j < g.fluid.NumY; j++ {
			if g.fluid.IsSolid(i, g.fluid.NumY-j-1) {
				idx := g.fluidToImageIndex(i, j)
				g.image.Pix[idx+0] = 0
				g.image.Pix[idx+1] = 0
				g.image.Pix[idx+2] = 0
				g.image.Pix[idx+3] = 0xff
			}
		}
	})

	// render
	g.eImage.ReplacePixels(g.image.Pix)
	screen.DrawImage(g.eImage, drawOpts)
	if g.showArrows {
		vel := g.fluid.Velocity()
		step := 10
		scale := float64(10)
		for i := 1; i < vel.NumX-1; i += step {
			for j := 1; j < vel.NumY-1; j += step {
				u, v, _ := vel.Value(i, j)
				x1 := float64(i)
				y1 := float64(g.fluid.NumY - j - 1)
				x2 := x1 + float64(u)*scale
				y2 := y1 - float64(v)*scale
				ebitenutil.DrawLine(screen, x1, y1, x2, y2, color.RGBA{0, 0, 0, 128})
			}
		}
	}
	show := "Pressure"
	if g.showSmoke {
		show = "Smoke"
	}
	wallInfo := fmt.Sprintf("U:%v D:%v L:%v R:%v",
		g.wallTop, g.wallBottom, g.wallLeft, g.wallRight)
	ebitenutil.DebugPrint(screen,
		fmt.Sprintf("FPS: %0.2f [S]how:%v Speed:%0.2f\n[J]et:%v [G]ravity:%0.2f [V]ort:%v [A]rrows:%v Paused:%v\n%s",
			ebiten.ActualFPS(), show, g.speed, g.jet, g.fluid.Gravity, g.fluid.Confinement != 0, g.showArrows, g.paused, wallInfo),
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

func (g *Game) applyWallSettings() {
	for i := 0; i < g.fluid.NumX; i++ {
		g.fluid.SetSolid(i, 0, g.wallBottom)
		g.fluid.SetSolid(i, g.fluid.NumY-1, g.wallTop)
	}
	for j := 0; j < g.fluid.NumY; j++ {
		g.fluid.SetSolid(0, j, g.wallLeft)
		g.fluid.SetSolid(g.fluid.NumX-1, j, g.wallRight)
	}
}

func (g *Game) clearWalls() {
	for i := 0; i < g.fluid.NumX; i++ {
		for j := 0; j < g.fluid.NumY; j++ {
			g.fluid.SetSolid(i, j, false)
		}
	}
	g.applyWallSettings()
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

package main

import (
	"fmt"
	"image"
	"image/color"
	"log"
	"math"
	"runtime"
	"sync"

	"github.com/TheFellow/fluid/pkg/fluid"
	"github.com/hajimehoshi/ebiten/v2"
	"github.com/hajimehoshi/ebiten/v2/ebitenutil"
	"github.com/hajimehoshi/ebiten/v2/inpututil"
	"github.com/hajimehoshi/ebiten/v2/vector"
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

// Tool modes for mouse interaction
type ToolMode int

const (
	ToolWall ToolMode = iota
	ToolErase
	ToolForce
	ToolSmoke
	ToolSource
	ToolSink
	toolCount // sentinel for cycling
)

func (t ToolMode) String() string {
	switch t {
	case ToolWall:
		return "Wall"
	case ToolErase:
		return "Erase"
	case ToolForce:
		return "Force"
	case ToolSmoke:
		return "Smoke"
	case ToolSource:
		return "Source"
	case ToolSink:
		return "Sink"
	}
	return "?"
}

// Visualization modes
type VizMode int

const (
	VizSmoke VizMode = iota
	VizPressure
	VizVelMag
	VizVorticity
	vizCount // sentinel for cycling
)

func (v VizMode) String() string {
	switch v {
	case VizSmoke:
		return "Smoke"
	case VizPressure:
		return "Pressure"
	case VizVelMag:
		return "Velocity"
	case VizVorticity:
		return "Vorticity"
	}
	return "?"
}

// Preset scenarios
type Preset int

const (
	PresetJet Preset = iota
	PresetCavity
	PresetKarman
)

func (p Preset) String() string {
	switch p {
	case PresetJet:
		return "Jet"
	case PresetCavity:
		return "Cavity"
	case PresetKarman:
		return "Karman"
	}
	return "?"
}

// Particle for Lagrangian tracer visualization
type Particle struct {
	X, Y   float32
	R, G, B uint8
	Age     float32
	MaxAge  float32
}

// Source region for persistent velocity injection
type Source struct {
	I, J int
	U, V float32
}

type Game struct {
	fluid  *fluid.Fluid
	image  *image.RGBA
	eImage *ebiten.Image

	// Visualization
	vizMode    VizMode
	showArrows bool
	showHelp   bool

	// Simulation control
	jet    bool
	paused bool
	speed  float32

	// Wall settings
	wallTop    bool
	wallBottom bool
	wallLeft   bool
	wallRight  bool

	// Tool state
	tool          ToolMode
	brushRadius   int
	brushStrength float32

	// Mouse tracking
	dragging  bool
	dragValue bool // for wall/erase tools
	prevMX    int  // previous mouse position (screen coords)
	prevMY    int
	prevI     int // previous fluid cell
	prevJ     int

	// Particles
	particles    []Particle
	showParticle bool
	hueAngle     float32 // rotating hue for injected particles
	maxParticles int

	// Sources and sinks
	sources []Source
	sinks   []Source

	// Preset
	preset Preset
}

func NewGame() *Game {
	f := fluid.New(1000, fluidWidth, fluidHeight, 1.0/100.0)
	img := image.NewRGBA(image.Rect(0, 0, fluidWidth+2, fluidHeight+2))
	g := &Game{
		fluid:         f,
		image:         img,
		eImage:        ebiten.NewImageFromImage(img),
		vizMode:       VizSmoke,
		speed:         1.0,
		wallTop:       true,
		wallBottom:    true,
		wallLeft:      true,
		wallRight:     false,
		tool:          ToolWall,
		brushRadius:   3,
		brushStrength: 1.0,
		showParticle:  true,
		maxParticles:  10000,
		preset:        PresetJet,
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
	g.handleKeyboard()
	g.handleMouse()
	g.applySources()

	if g.jet {
		for j := fluidHeight/2 - 100; j < fluidHeight/2+100; j++ {
			g.fluid.SetVelocity(1, j, 4.0, 0)
			g.fluid.AddSmoke(1, j, 1.0)
		}
	}

	dt := float32(1.0/120.0) * g.speed
	if !g.paused {
		g.fluid.Simulate(dt)
		g.advectParticles(dt)
	}
	return nil
}

func (g *Game) handleKeyboard() {
	// Visualization
	if inpututil.IsKeyJustPressed(ebiten.KeyS) {
		g.vizMode = (g.vizMode + 1) % vizCount
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyA) {
		g.showArrows = !g.showArrows
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyH) {
		g.showHelp = !g.showHelp
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyP) {
		g.showParticle = !g.showParticle
	}

	// Simulation
	if inpututil.IsKeyJustPressed(ebiten.KeyJ) {
		g.jet = !g.jet
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyV) {
		if g.fluid.Confinement != 0 {
			g.fluid.Confinement = 0
		} else {
			g.fluid.Confinement = 0.1
		}
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyB) {
		g.fluid.UseBFECC = !g.fluid.UseBFECC
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
		g.resetToPreset()
	}

	// Walls
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

	// Tool selection
	if inpututil.IsKeyJustPressed(ebiten.KeyW) {
		g.tool = ToolWall
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyE) {
		g.tool = ToolErase
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyF) {
		g.tool = ToolForce
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyK) {
		g.tool = ToolSmoke
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyO) {
		g.tool = ToolSource
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyN) {
		g.tool = ToolSink
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyTab) {
		g.tool = (g.tool + 1) % toolCount
	}

	// Brush adjustments
	if inpututil.IsKeyJustPressed(ebiten.KeyBracketRight) {
		g.brushRadius++
		if g.brushRadius > 20 {
			g.brushRadius = 20
		}
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyBracketLeft) {
		g.brushRadius--
		if g.brushRadius < 1 {
			g.brushRadius = 1
		}
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyPeriod) {
		g.brushStrength *= 1.5
		if g.brushStrength > 50 {
			g.brushStrength = 50
		}
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyComma) {
		g.brushStrength /= 1.5
		if g.brushStrength < 0.1 {
			g.brushStrength = 0.1
		}
	}

	// Presets
	if inpututil.IsKeyJustPressed(ebiten.KeyDigit1) {
		g.preset = PresetJet
		g.resetToPreset()
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyDigit2) {
		g.preset = PresetCavity
		g.resetToPreset()
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyDigit3) {
		g.preset = PresetKarman
		g.resetToPreset()
	}
}

func (g *Game) handleMouse() {
	mx, my := ebiten.CursorPosition()

	// Right-drag: always applies force impulse
	if ebiten.IsMouseButtonPressed(ebiten.MouseButtonRight) {
		i, j := g.screenToFluid(mx, my)
		if g.inBounds(i, j) {
			// Compute velocity from mouse delta
			dx := float32(mx-g.prevMX) * g.brushStrength * 0.5
			dy := -float32(my-g.prevMY) * g.brushStrength * 0.5 // flip Y
			if dx != 0 || dy != 0 {
				g.fluid.ApplyForceRadius(i, j, dx, dy, g.brushRadius)
			}
		}
	}

	// Left-drag: applies active tool
	if inpututil.IsMouseButtonJustPressed(ebiten.MouseButtonLeft) {
		i, j := g.screenToFluid(mx, my)
		if g.inBounds(i, j) {
			g.dragging = true
			switch g.tool {
			case ToolWall:
				g.dragValue = !g.fluid.IsSolid(i, j)
				g.drawLine(i, j, i, j, g.dragValue)
			case ToolErase:
				g.dragValue = false
				g.drawLine(i, j, i, j, false)
			case ToolForce:
				// nothing on click, handled during drag
			case ToolSmoke:
				g.applySmokeBrush(i, j)
			case ToolSource:
				g.sources = append(g.sources, Source{I: i, J: j, U: 5.0 * g.brushStrength, V: 0})
			case ToolSink:
				g.sinks = append(g.sinks, Source{I: i, J: j, U: -5.0 * g.brushStrength, V: 0})
			}
			g.prevI = i
			g.prevJ = j
		}
	} else if ebiten.IsMouseButtonPressed(ebiten.MouseButtonLeft) && g.dragging {
		i, j := g.screenToFluid(mx, my)
		if g.inBounds(i, j) {
			switch g.tool {
			case ToolWall:
				if i != g.prevI || j != g.prevJ {
					g.drawLine(g.prevI, g.prevJ, i, j, g.dragValue)
				}
			case ToolErase:
				if i != g.prevI || j != g.prevJ {
					g.drawLine(g.prevI, g.prevJ, i, j, false)
				}
			case ToolForce:
				dx := float32(mx-g.prevMX) * g.brushStrength * 0.5
				dy := -float32(my-g.prevMY) * g.brushStrength * 0.5
				if dx != 0 || dy != 0 {
					g.fluid.ApplyForceRadius(i, j, dx, dy, g.brushRadius)
				}
			case ToolSmoke:
				g.applySmokeBrush(i, j)
			}
			g.prevI = i
			g.prevJ = j
		}
	} else if g.dragging {
		g.dragging = false
	}

	g.prevMX = mx
	g.prevMY = my
}

func (g *Game) applySmokeBrush(ci, cj int) {
	r := g.brushRadius
	for i := ci - r; i <= ci+r; i++ {
		for j := cj - r; j <= cj+r; j++ {
			if !g.inBounds(i, j) {
				continue
			}
			dx := float32(i - ci)
			dy := float32(j - cj)
			if dx*dx+dy*dy > float32(r*r) {
				continue
			}
			g.fluid.AddSmoke(i, j, g.brushStrength*0.3)
		}
	}
	// Spawn particles at cursor
	if g.showParticle && len(g.particles) < g.maxParticles {
		g.spawnParticles(ci, cj, 5)
	}
}

func (g *Game) applySources() {
	for _, src := range g.sources {
		if g.inBounds(src.I, src.J) && !g.fluid.IsSolid(src.I, src.J) {
			g.fluid.SetVelocity(src.I, src.J, src.U, src.V)
			g.fluid.AddSmoke(src.I, src.J, 0.5)
		}
	}
	for _, sink := range g.sinks {
		if g.inBounds(sink.I, sink.J) && !g.fluid.IsSolid(sink.I, sink.J) {
			g.fluid.SetVelocity(sink.I, sink.J, sink.U, sink.V)
		}
	}
}

// Particle management

func (g *Game) spawnParticles(ci, cj, count int) {
	h := g.fluid.H()
	for k := 0; k < count; k++ {
		if len(g.particles) >= g.maxParticles {
			break
		}
		// World position at cell center
		x := (float32(ci) + 0.5) * h
		y := (float32(cj) + 0.5) * h
		r, gr, b := hueToRGB(g.hueAngle)
		g.particles = append(g.particles, Particle{
			X: x, Y: y,
			R: r, G: gr, B: b,
			Age: 0, MaxAge: 15.0,
		})
	}
	g.hueAngle += 2.0
	if g.hueAngle >= 360 {
		g.hueAngle -= 360
	}
}

func (g *Game) advectParticles(dt float32) {
	if !g.showParticle {
		return
	}
	h := g.fluid.H()
	alive := g.particles[:0]
	for i := range g.particles {
		p := &g.particles[i]
		p.Age += dt

		if p.Age > p.MaxAge {
			continue
		}

		// RK2 midpoint integration
		u1, v1 := g.fluid.SampleVelocity(p.X, p.Y)
		midX := p.X + 0.5*dt*u1
		midY := p.Y + 0.5*dt*v1
		u2, v2 := g.fluid.SampleVelocity(midX, midY)
		p.X += dt * u2
		p.Y += dt * v2

		// Check bounds
		fi := int(p.X / h)
		fj := int(p.Y / h)
		if fi < 0 || fi >= g.fluid.NumX || fj < 0 || fj >= g.fluid.NumY {
			continue
		}
		if g.fluid.IsSolid(fi, fj) {
			continue
		}
		alive = append(alive, *p)
	}
	g.particles = alive
}

var drawOpts = &ebiten.DrawImageOptions{Filter: ebiten.FilterLinear}

func (g *Game) Draw(screen *ebiten.Image) {
	// Render the selected scalar field
	switch g.vizMode {
	case VizSmoke:
		g.drawScalarField(g.fluid.Smoke())
	case VizPressure:
		g.drawScalarField(g.fluid.Pressure())
	case VizVelMag:
		g.drawScalarField(g.fluid.VelocityMagnitude())
	case VizVorticity:
		g.drawVorticityField(g.fluid.Vorticity())
	}

	// Overlay solid cells as black
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

	// Draw source/sink indicators
	for _, src := range g.sources {
		sj := g.fluid.NumY - src.J - 1
		idx := g.fluidToImageIndex(src.I, sj)
		if idx >= 0 && idx+3 < len(g.image.Pix) {
			g.image.Pix[idx+0] = 0
			g.image.Pix[idx+1] = 255
			g.image.Pix[idx+2] = 0
			g.image.Pix[idx+3] = 255
		}
	}
	for _, sink := range g.sinks {
		sj := g.fluid.NumY - sink.J - 1
		idx := g.fluidToImageIndex(sink.I, sj)
		if idx >= 0 && idx+3 < len(g.image.Pix) {
			g.image.Pix[idx+0] = 255
			g.image.Pix[idx+1] = 0
			g.image.Pix[idx+2] = 0
			g.image.Pix[idx+3] = 255
		}
	}

	g.eImage.WritePixels(g.image.Pix)
	screen.DrawImage(g.eImage, drawOpts)

	// Velocity arrows / streamlines
	if g.showArrows {
		g.drawVelocityOverlay(screen)
	}

	// Particles
	if g.showParticle {
		g.drawParticles(screen)
	}

	// HUD
	g.drawHUD(screen)

	// Help overlay
	if g.showHelp {
		g.drawHelp(screen)
	}
}

func (g *Game) drawScalarField(sf fluid.ScalarField) {
	parallelRange(0, sf.NumX, func(i int) {
		for j := 0; j < sf.NumY; j++ {
			p, err := sf.Value(i, sf.NumY-j-1)
			if err != nil {
				continue
			}
			c := getSciValue(p, sf.MinValue, sf.MaxValue)
			idx := g.fluidToImageIndex(i, j)
			g.image.Pix[idx+0] = c.R
			g.image.Pix[idx+1] = c.G
			g.image.Pix[idx+2] = c.B
			g.image.Pix[idx+3] = c.A
		}
	})
}

func (g *Game) drawVorticityField(sf fluid.ScalarField) {
	parallelRange(0, sf.NumX, func(i int) {
		for j := 0; j < sf.NumY; j++ {
			v, err := sf.Value(i, sf.NumY-j-1)
			if err != nil {
				continue
			}
			c := getDivergingColor(v, sf.MinValue, sf.MaxValue)
			idx := g.fluidToImageIndex(i, j)
			g.image.Pix[idx+0] = c.R
			g.image.Pix[idx+1] = c.G
			g.image.Pix[idx+2] = c.B
			g.image.Pix[idx+3] = c.A
		}
	})
}

func (g *Game) drawVelocityOverlay(screen *ebiten.Image) {
	vel := g.fluid.Velocity()
	step := 10
	scale := float32(10)
	for i := 1; i < vel.NumX-1; i += step {
		for j := 1; j < vel.NumY-1; j += step {
			u, v, _ := vel.Value(i, j)
			x1 := float32(i)
			y1 := float32(g.fluid.NumY - j - 1)
			x2 := x1 + u*scale
			y2 := y1 - v*scale
			vector.StrokeLine(screen, x1, y1, x2, y2, 1, color.RGBA{G: 200, A: 160}, false)
		}
	}
}

func (g *Game) drawParticles(screen *ebiten.Image) {
	h := g.fluid.H()
	for i := range g.particles {
		p := &g.particles[i]
		sx := p.X / h
		sy := float32(g.fluid.NumY) - p.Y/h
		// Fade out as particle ages
		alpha := uint8(255 * (1.0 - p.Age/p.MaxAge))
		if alpha < 20 {
			alpha = 20
		}
		vector.DrawFilledRect(screen, sx-1, sy-1, 2, 2, color.RGBA{R: p.R, G: p.G, B: p.B, A: alpha}, false)
	}
}

func (g *Game) drawHUD(screen *ebiten.Image) {
	bfecc := "Off"
	if g.fluid.UseBFECC {
		bfecc = "On"
	}
	vort := "Off"
	if g.fluid.Confinement != 0 {
		vort = "On"
	}
	wallInfo := fmt.Sprintf("U:%v D:%v L:%v R:%v",
		g.wallTop, g.wallBottom, g.wallLeft, g.wallRight)
	ebitenutil.DebugPrint(screen,
		fmt.Sprintf("FPS: %.0f  [S]:%s  Speed:%.2f  [B]FECC:%s\n"+
			"[J]et:%v  [V]ort:%s  [A]rrows:%v  [P]articles:%v\n"+
			"Tool[Tab]:%s  Brush:[/] %d  Str:<,> %.1f\n"+
			"Preset[1-3]:%s  %s  [H]elp",
			ebiten.ActualFPS(), g.vizMode, g.speed, bfecc,
			g.jet, vort, g.showArrows, g.showParticle,
			g.tool, g.brushRadius, g.brushStrength,
			g.preset, wallInfo),
	)
}

func (g *Game) drawHelp(screen *ebiten.Image) {
	help := `
  === Controls ===
  S         Cycle viz: Smoke/Pressure/Velocity/Vorticity
  J         Toggle jet
  V         Toggle vorticity confinement
  B         Toggle BFECC advection
  A         Toggle velocity arrows
  P         Toggle particles
  Space     Pause/resume
  +/-       Speed up/down
  R         Reset to current preset
  Arrows    Toggle border walls
  C         Clear all interior walls

  === Tools (Tab to cycle) ===
  W         Wall tool
  E         Erase tool
  F         Force tool
  K         Smoke tool
  O         Source tool
  N         Sink tool
  [ / ]     Brush radius
  , / .     Brush strength

  === Mouse ===
  Left-drag   Apply active tool
  Right-drag  Force impulse (always)

  === Presets ===
  1         Jet flow
  2         Lid-driven cavity
  3         Karman vortex street

  H         Close this help
`
	ebitenutil.DebugPrintAt(screen, help, 10, 60)
}

func (g *Game) resetToPreset() {
	g.fluid.Reset()
	g.particles = g.particles[:0]
	g.sources = g.sources[:0]
	g.sinks = g.sinks[:0]

	// Clear all walls first
	for i := 0; i < g.fluid.NumX; i++ {
		for j := 0; j < g.fluid.NumY; j++ {
			g.fluid.SetSolid(i, j, false)
		}
	}

	switch g.preset {
	case PresetJet:
		g.wallTop = true
		g.wallBottom = true
		g.wallLeft = true
		g.wallRight = false
		g.jet = true

	case PresetCavity:
		g.wallTop = true
		g.wallBottom = true
		g.wallLeft = true
		g.wallRight = true
		g.jet = false
		// Top wall drives flow rightward via source row
		for i := 2; i < g.fluid.NumX-2; i++ {
			g.sources = append(g.sources, Source{
				I: i, J: g.fluid.NumY - 2,
				U: 3.0, V: 0,
			})
		}

	case PresetKarman:
		g.wallTop = true
		g.wallBottom = true
		g.wallLeft = true
		g.wallRight = false
		g.jet = true
		// Circular obstacle
		cx := g.fluid.NumX / 4
		cy := g.fluid.NumY / 2
		g.fluid.SetCircularObstacle(cx, cy, 8)
	}

	g.applyWallSettings()
}

func (g *Game) fluidToImageIndex(i, j int) int {
	return 4*i + j*g.image.Stride
}

func (g *Game) screenToFluid(mx, my int) (int, int) {
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

// hueToRGB converts a hue angle (0-360) to RGB.
func hueToRGB(h float32) (uint8, uint8, uint8) {
	h = float32(math.Mod(float64(h), 360))
	c := float32(1.0)
	x := c * (1.0 - float32(math.Abs(math.Mod(float64(h/60), 2)-1)))
	var r, g, b float32
	switch {
	case h < 60:
		r, g, b = c, x, 0
	case h < 120:
		r, g, b = x, c, 0
	case h < 180:
		r, g, b = 0, c, x
	case h < 240:
		r, g, b = 0, x, c
	case h < 300:
		r, g, b = x, 0, c
	default:
		r, g, b = c, 0, x
	}
	return uint8(r * 255), uint8(g * 255), uint8(b * 255)
}

func main() {
	ebiten.SetWindowSize(screenWidth, screenHeight)
	ebiten.SetWindowTitle("FluidSim")

	if err := ebiten.RunGame(NewGame()); err != nil {
		log.Fatal(err)
	}
}

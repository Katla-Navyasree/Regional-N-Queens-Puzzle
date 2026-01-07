import tkinter as tk
import random
import colorsys

# ---------- Configuration ----------
MIN_CELL_SIZE = 50    # minimum size for each square
MAX_CANVAS_SIZE = 600 # max canvas side in pixels

# ---------- Utilities ----------

def generate_palette(n):
    out = []
    for k in range(n):
        h = k / max(1, n)
        s = 0.55
        v = 0.95
        r, g, b = colorsys.hsv_to_rgb(h, s, v)
        out.append('#{:02x}{:02x}{:02x}'.format(int(r*255), int(g*255), int(b*255)))
    return out


def generate_blob_regions(n, max_restarts=200, min_seed_spacing=1):
    """
    Generate N contiguous irregular blobs (region ids 0..n-1) covering an n x n grid.
    Multi-source randomized growth ensures compact connected regions.
    """
    dirs = [(1,0),(-1,0),(0,1),(0,-1)]
    for attempt in range(max_restarts):
        region = [[-1]*n for _ in range(n)]
        seeds = []
        tries = 0
        while len(seeds) < n and tries < n*n*10:
            r = random.randrange(n)
            c = random.randrange(n)
            if region[r][c] != -1:
                tries += 1
                continue
            ok = True
            for (sr, sc) in seeds:
                if abs(sr - r) + abs(sc - c) < min_seed_spacing:
                    ok = False
                    break
            if ok:
                idx = len(seeds)
                region[r][c] = idx
                seeds.append((r, c))
            tries += 1
        for i in range(n):
            for j in range(n):
                if len(seeds) >= n:
                    break
                if region[i][j] == -1:
                    idx = len(seeds)
                    region[i][j] = idx
                    seeds.append((i,j))
            if len(seeds) >= n:
                break

        frontiers = {k: [seeds[k]] for k in range(n)}
        unfilled = n*n - n
        while unfilled > 0 and any(frontiers[k] for k in frontiers):
            rid = random.choice([k for k in frontiers if frontiers[k]])
            ci = random.randrange(len(frontiers[rid]))
            cr, cc = frontiers[rid].pop(ci)
            random.shuffle(dirs)
            for di, dj in dirs:
                nr, nc = cr + di, cc + dj
                if 0 <= nr < n and 0 <= nc < n and region[nr][nc] == -1:
                    region[nr][nc] = rid
                    frontiers[rid].append((nr, nc))
                    unfilled -= 1
                    break
        if unfilled == 0:
            return region
    region = [[j % n for j in range(n)] for i in range(n)]
    return region

# ---------- Conflict detection ----------

def detect_conflicts(queens, region):
    n = len(region)
    rowmap = {}
    colmap = {}
    d1map = {}
    d2map = {}
    regionmap = {}
    for (r,c) in queens:
        rowmap.setdefault(r, []).append((r,c))
        colmap.setdefault(c, []).append((r,c))
        d1map.setdefault(r-c, []).append((r,c))
        d2map.setdefault(r+c, []).append((r,c))
        rid = region[r][c]
        regionmap.setdefault(rid, []).append((r,c))

    conflicts_positions = set()
    row_conflicts = set(k for k,lst in rowmap.items() if len(lst) > 1)
    for r in row_conflicts:
        conflicts_positions.update(rowmap[r])
    col_conflicts = set(k for k,lst in colmap.items() if len(lst) > 1)
    for c in col_conflicts:
        conflicts_positions.update(colmap[c])
    d1_conflicts = set(k for k,lst in d1map.items() if len(lst) > 1)
    for k in d1_conflicts:
        conflicts_positions.update(d1map[k])
    d2_conflicts = set(k for k,lst in d2map.items() if len(lst) > 1)
    for k in d2_conflicts:
        conflicts_positions.update(d2map[k])
    region_conflicted_ids = set(k for k,lst in regionmap.items() if len(lst) > 1)

    return conflicts_positions, region_conflicted_ids, row_conflicts, col_conflicts, d1_conflicts, d2_conflicts

# ---------- Line-rectangle clipping (Liang–Barsky) ----------

def liang_barsky_clip(x0, y0, x1, y1, xmin, ymin, xmax, ymax):
    dx = x1 - x0
    dy = y1 - y0
    p = [-dx, dx, -dy, dy]
    q = [x0 - xmin, xmax - x0, y0 - ymin, ymax - y0]
    u1 = 0.0
    u2 = 1.0
    for pi, qi in zip(p, q):
        if pi == 0:
            if qi < 0:
                return None
        else:
            t = qi / pi
            if pi < 0:
                if t > u2:
                    return None
                if t > u1:
                    u1 = t
            else:
                if t < u1:
                    return None
                if t < u2:
                    u2 = t
    if u1 > u2:
        return None
    cx0 = x0 + u1 * dx
    cy0 = y0 + u1 * dy
    cx1 = x0 + u2 * dx
    cy1 = y0 + u2 * dy
    return cx0, cy0, cx1, cy1

# ---------- Solver ----------

def solve_with_colors(region):
    n = len(region)
    board = [[0]*n for _ in range(n)]
    cols = [False]*n
    diag1 = [False]*(2*n-1)
    diag2 = [False]*(2*n-1)
    colors_used = [False]*n

    solution = None

    def backtrack(r):
        nonlocal solution
        if r == n:
            solution = [row[:] for row in board]
            return True
        cols_order = list(range(n))
        random.shuffle(cols_order)
        for c in cols_order:
            color = region[r][c]
            d1 = r - c + (n-1)
            d2 = r + c
            if not cols[c] and not diag1[d1] and not diag2[d2] and not colors_used[color]:
                board[r][c] = 1
                cols[c] = diag1[d1] = diag2[d2] = colors_used[color] = True
                if backtrack(r+1):
                    return True
                board[r][c] = 0
                cols[c] = diag1[d1] = diag2[d2] = False
                colors_used[color] = False
        return False

    if backtrack(0):
        return solution
    return None

# ---------- GUI App ----------

class QueensApp:
    def __init__(self, master, n=8):
        self.master = master
        self.n = n
        self.update_canvas_size()

        self.canvas = tk.Canvas(master, width=self.canvas_size, height=self.canvas_size, bg='white')
        self.canvas.grid(row=0, column=0, columnspan=6)

        # --- Size selector ---
        tk.Label(master, text="Board Size:").grid(row=1, column=0)
        self.size_var = tk.IntVar(value=n)
        size_menu = tk.OptionMenu(master, self.size_var, *list(range(4, 13)), command=self.change_size)
        size_menu.grid(row=1, column=1)

        # --- Buttons ---
        self.btn_check = tk.Button(master, text='Check Solution', command=self.check_solution)
        self.btn_check.grid(row=1, column=2)
        self.btn_auto = tk.Button(master, text='Auto Solve', command=self.auto_solve)
        self.btn_auto.grid(row=1, column=3)
        self.btn_reset = tk.Button(master, text='Reset', command=self.reset_board)
        self.btn_reset.grid(row=1, column=4)
        self.btn_new = tk.Button(master, text='New Regions', command=self.new_regions)
        self.btn_new.grid(row=1, column=5)

        self.status = tk.Label(master, text='Place queens (♛) and watch conflicts highlight.')
        self.status.grid(row=2, column=0, columnspan=6)

        self.generate_solvable_region()
        self.queens = []

        self.canvas.bind('<Button-1>', self.on_click)
        self.draw()

    def update_canvas_size(self):
        """Recompute canvas + cell size dynamically."""
        self.cell_size = max(MIN_CELL_SIZE, MAX_CANVAS_SIZE // self.n)
        self.canvas_size = self.n * self.cell_size

    def generate_solvable_region(self):
        tries = 0
        while True:
            self.region = generate_blob_regions(self.n)
            if solve_with_colors(self.region) is not None:
                break
            tries += 1
            if tries > 200:
                break
        self.palette = generate_palette(self.n)

    def change_size(self, new_size):
        self.n = int(new_size)
        self.update_canvas_size()
        self.canvas.config(width=self.canvas_size, height=self.canvas_size)
        self.generate_solvable_region()
        self.queens = []
        self.status.config(text=f'New {self.n}x{self.n} solvable board created.')
        self.draw()

    def draw(self):
        self.canvas.delete('all')
        for r in range(self.n):
            for c in range(self.n):
                rid = self.region[r][c]
                color = self.palette[rid % len(self.palette)]
                x1, y1 = c*self.cell_size, r*self.cell_size
                x2, y2 = x1 + self.cell_size, y1 + self.cell_size
                self.canvas.create_rectangle(x1, y1, x2, y2, fill=color, outline='black')
        self.highlight_and_draw()

    def draw_stripes_over_cells(self, cells):
        if not cells:
            return
        n = self.n
        cs = self.cell_size
        min_r = min(r for r, _ in cells)
        max_r = max(r for r, _ in cells)
        min_c = min(c for _, c in cells)
        max_c = max(c for _, c in cells)
        bx1 = min_c * cs
        by1 = min_r * cs
        bx2 = (max_c + 1) * cs
        by2 = (max_r + 1) * cs
        w = bx2 - bx1
        step = max(6, int(cs/6))

        rects = []
        for (r, c) in cells:
            cx1 = c * cs
            cy1 = r * cs
            cx2 = cx1 + cs
            cy2 = cy1 + cs
            rects.append((cx1, cy1, cx2, cy2))

        for offset in range(-w, w*2, step):
            x_start = bx1 + offset
            y_start = by1
            x_end = x_start - w
            y_end = by2
            for (cx1, cy1, cx2, cy2) in rects:
                seg = liang_barsky_clip(x_start, y_start, x_end, y_end, cx1, cy1, cx2, cy2)
                if seg:
                    sx0, sy0, sx1, sy1 = seg
                    self.canvas.create_line(sx0, sy0, sx1, sy1, fill='#ff4d4d', width=3)

    def highlight_and_draw(self):
        conflicts_positions, region_conflicted_ids, row_conflicts, col_conflicts, d1_conflicts, d2_conflicts = detect_conflicts(self.queens, self.region)
        for rid in region_conflicted_ids:
            cells = [(r, c) for r in range(self.n) for c in range(self.n) if self.region[r][c] == rid]
            self.draw_stripes_over_cells(cells)
        for r in row_conflicts:
            cells = [(r, c) for c in range(self.n)]
            self.draw_stripes_over_cells(cells)
        for c in col_conflicts:
            cells = [(r, c) for r in range(self.n)]
            self.draw_stripes_over_cells(cells)
        for k in d1_conflicts:
            cells = [(r, r - k) for r in range(self.n) if 0 <= (r - k) < self.n]
            self.draw_stripes_over_cells(cells)
        for k in d2_conflicts:
            cells = [(r, k - r) for r in range(self.n) if 0 <= (k - r) < self.n]
            self.draw_stripes_over_cells(cells)
        for (r, c) in self.queens:
            x, y = c * self.cell_size + self.cell_size // 2, r * self.cell_size + self.cell_size // 2
            self.canvas.create_text(x, y, text='♛', font=('Arial', max(12, int(self.cell_size * 0.6))), fill='black')

    def on_click(self, event):
        c = event.x // self.cell_size
        r = event.y // self.cell_size
        if not (0 <= r < self.n and 0 <= c < self.n):
            return
        if (r,c) in self.queens:
            self.queens.remove((r,c))
            self.status.config(text='Removed queen.')
        else:
            self.queens.append((r,c))
            self.status.config(text='Placed queen.')
        self.draw()

    def check_solution(self):
        conflicts_positions, region_conflicted_ids, *_ = detect_conflicts(self.queens, self.region)
        if conflicts_positions or region_conflicted_ids:
            self.status.config(text='❌ Conflicts detected — highlighted on board.')
        elif len(self.queens) != self.n:
            self.status.config(text=f'❌ Place exactly {self.n} queens (currently {len(self.queens)}).')
        else:
            self.status.config(text='🎉 You solved it!')
        self.draw()

    def reset_board(self):
        self.queens = []
        self.status.config(text='Board reset. Place queens again.')
        self.draw()

    def new_regions(self):
        self.generate_solvable_region()
        self.queens = []
        self.status.config(text=f'New regions generated for {self.n}x{self.n} board.')
        self.draw()

    def auto_solve(self):
        sol = solve_with_colors(self.region)
        if sol is None:
            self.status.config(text='No solution found for this partition.')
            return
        self.queens = []
        for r in range(self.n):
            for c in range(self.n):
                if sol[r][c] == 1:
                    self.queens.append((r,c))
        self.status.config(text='🤖 Auto-solved!')
        self.draw()


if __name__ == '__main__':
    root = tk.Tk()
    root.title("N-Queens with Regions")
    app = QueensApp(root, n=7)
    root.mainloop()

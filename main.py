#!/usr/bin/env python3
import sys
import tkinter as tk
import math
from pointing import az_el_from_geodetic  # device-agnostic pointing math

def runPredictionTool(checkbox_dict, tle_dict, my_lat, my_lon):
    """
    Called when you press “Run Prediction.” Gathers the checked satellites,
    then builds the 2×2 figure and starts the animation.
    """
    # 1) Gather satellites with IntVar == 1
    selected = {name: tle_dict[name] for name, var in checkbox_dict.items() if var.get() == 1}
    if not selected:
        from tkinter import messagebox
        messagebox.showwarning("No Satellites Selected", "Please check at least one satellite before running.")
        return

    # 2) Imports used only by the plotting/animation
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.patheffects as pe

    from mpl_toolkits.basemap import Basemap
    from matplotlib import animation
    from datetime import datetime, timedelta
    from collections import deque

    from tle_parser import parse_and_convert_tle
    from kep_to_state import ConvertKepToStateVectors

    # Build time window for your TLE conversion
    now   = datetime.utcnow()
    start = now.strftime('%Y %m %d %H %M %S')
    end   = (now + timedelta(minutes=90)).strftime('%Y %m %d %H %M %S')

    # TLE → Kepler params and time vector (keeps your existing pipeline)
    kepler_dict, time_vec, epoch = parse_and_convert_tle("amateur.tle", start, end)

    # ────────────────────────────────────────────────────────────────────────
    # Figure + layout (2×2) — bottom-left is split into two polar subplots
    fig = plt.figure(figsize=(14, 7))
    gs = fig.add_gridspec(2, 2, hspace=0.25, wspace=0.25)

    # Top row: maps
    ax1 = fig.add_subplot(gs[0, 0])   # Global
    ax2 = fig.add_subplot(gs[0, 1])   # Near-sided

    # Bottom-left split into two POLAR plots (Azimuth + Elevation)
    subgs = gs[1, 0].subgridspec(2, 1, hspace=0.35)
    ax_az = fig.add_subplot(subgs[0, 0], projection='polar')  # Compass rose
    ax_el = fig.add_subplot(subgs[1, 0], projection='polar')  # Elevation semicircle

    # Bottom-right: altitude panel
    ax4 = fig.add_subplot(gs[1, 1])

    # Compass label helper (16-point compass)
    def az_to_compass(az):
        dirs = ['N','NNE','NE','ENE','E','ESE','SE','SSE',
                'S','SSW','SW','WSW','W','WNW','NW','NNW']
        return dirs[int((az/22.5)+0.5) % 16]

    # ── Style helpers for the two polar gauges
    def init_az_compass(ax):
        ax.set_facecolor('black')
        ax.set_theta_zero_location('N')  # 0° at North
        ax.set_theta_direction(-1)       # clockwise azimuth
        ax.set_rlim(0, 1.0)
        ax.set_rticks([])
        ax.set_xticklabels([])
        # ax.set_title("Azimuth", color='red', pad=8)
        # Title ABOVE the compass (outlined like Elevation)
        ax.text(0.5, 1.08, "Azimuth",
                transform=ax.transAxes, ha='center', va='bottom',
                color='white', fontsize=12, clip_on=False,
                path_effects=[pe.withStroke(linewidth=3, foreground='black')])

        # Compass rings and spokes
        for r in (0.33, 0.66, 1.0):
            ax.plot([0, 2*math.pi], [r, r], color='white', alpha=0.15, linewidth=1)
        for ang in range(0, 360, 30):
            t = math.radians(ang)
            ax.plot([t, t], [0.0, 1.0], color='white', alpha=0.15, linewidth=1)
        # NSEW labels (slightly inside so they don't clip)
        for ang, lab in [(0, 'N'), (90, 'E'), (180, 'S'), (270, 'W')]:
            ax.text(math.radians(ang), 0.7, lab, color='red',
                    ha='center', va='bottom', fontsize=10)

    def init_el_gauge(ax):
        ax.set_facecolor('black')
        ax.set_theta_zero_location('W')  # 0° at left
        ax.set_theta_direction(-1)       # clockwise toward up (90° at top)
        ax.set_thetamin(0); ax.set_thetamax(90)  # show only 0..90°
        ax.set_rlim(0, 1.0)
        ax.set_rticks([])
        ax.set_xticklabels([])
        ax_el.text(0.5, -0.18, "Elevation",
           transform=ax_el.transAxes, ha='center', va='top',
           color='white', fontsize=12, clip_on=False,
           path_effects=[pe.withStroke(linewidth=3, foreground='black')])

        # Subtle arc guides
        for r in (0.55, 0.8, 1.0):
            ax.plot([math.radians(0), math.radians(90)], [r, r],
                    color='white', alpha=0.18, linewidth=1)

        # Major ticks and labels (high contrast + outline)
        major_degs = (0, 10, 30, 60, 90)
        for ang in major_degs:
            t = math.radians(ang)
            # tick mark
            ax.plot([t, t], [0.94, 1.0], color='white', linewidth=1.6, alpha=0.9, zorder=4)
            # label slightly outside the rim; outline makes it readable over anything
            ax.text(t, 1.06, f"{ang}°",
                    color='white', fontsize=12, fontweight='bold',
                    ha='center', va='bottom', clip_on=False,
                    path_effects=[pe.withStroke(linewidth=3, foreground='black')])

        # lightly shade 0–10° mask
        ax.bar(math.radians(5), 1.0, width=math.radians(10), bottom=0.0,
            alpha=0.14, color='red', edgecolor=None)

    # def init_el_gauge(ax):
    #     ax.set_facecolor('black')
    #     ax.set_theta_zero_location('W')  # 0° at left
    #     ax.set_theta_direction(-1)       # clockwise toward up (90° at top)
    #     ax.set_thetamin(0); ax.set_thetamax(90)  # show only 0..90°
    #     ax.set_rlim(0, 1.0)
    #     ax.set_rticks([])
    #     ax.set_xticklabels([])
    #     ax.set_title("Elevation", color='red', pad=8)
    #     # arc guides + ticks
    #     for r in (0.55, 0.8, 1.0):
    #         ax.plot([math.radians(0), math.radians(90)], [r, r], color='white', alpha=0.15, linewidth=1)
    #     for ang in (0, 10, 30, 60, 90):
    #         t = math.radians(ang)
    #         ax.plot([t, t], [0.0, 1.0], color='white', alpha=0.15, linewidth=1)
    #         ax.text(t, 0.98, f"{ang}°", color='white', ha='center', va='bottom', fontsize=9)
    #     # lightly shade 0–10° mask
    #     ax.bar(math.radians(5), 1.0, width=math.radians(10), bottom=0.0, alpha=0.12, color='red', edgecolor=None)

    # ── Top-Left: Global View
    map1 = Basemap(
        projection='mill',
        llcrnrlat=-90, urcrnrlat=90,
        llcrnrlon=-180, urcrnrlon=180,
        resolution='c',
        ax=ax1
    )
    map1.drawmapboundary(fill_color='black')
    map1.fillcontinents(color='gray', lake_color='black')
    map1.drawcoastlines()
    ax1.set_facecolor('black')
    ax1.set_title("Global View", color='white')

    # Plot QTH once on global
    x_qth, y_qth = map1(my_lon, my_lat)
    map1.plot(x_qth, y_qth, 'go', markersize=8)
    ax1.annotate('Me', xy=(x_qth, y_qth), xytext=(x_qth + 5, y_qth + 5), color='white')

    # ── Top-Right: Near-Sided Perspective
    initial_height = 2000  # km (will update each frame)
    map2 = Basemap(
        projection='nsper',
        lon_0=my_lon,
        lat_0=my_lat,
        satellite_height=initial_height * 1000.0,
        resolution='l',
        ax=ax2
    )
    map2.drawcoastlines()
    map2.fillcontinents(color='coral', lake_color='aqua')
    map2.drawparallels(np.arange(-90, 120, 30))
    map2.drawmeridians(np.arange(0, 420, 60))
    map2.drawmapboundary(fill_color='aqua')
    ax2.set_facecolor('black')
    ax2.set_title(f"Near-Sided Perspective @ {initial_height} km", color='white')

    # ── Bottom-Right: Altitude vs. Time
    ax4.set_facecolor('darkslategray')
    ax4.set_title("Altitude vs. Time (rolling)", color='white')
    ax4.set_xlabel("Frame #", color='white')
    ax4.set_ylabel("Altitude (km)", color='white')
    ax4.tick_params(colors='white')
    max_history = 100
    alt_hist = deque(maxlen=max_history)
    frame_hist = deque(maxlen=max_history)
    line4, = ax4.plot([], [], 'y-', linewidth=2)

    # ────────────────────────────────────────────────────────────────────────
    # Utility to split “wrapping” ground tracks on map edges
    def split_segments(lons, lats, max_jump=180.0):
        segs = []
        start = 0
        for i in range(1, len(lons)):
            if abs(lons[i] - lons[i - 1]) > max_jump:
                segs.append((lons[start:i], lats[start:i]))
                start = i
        segs.append((lons[start:], lats[start:]))
        return segs

    track_objs = []  # Matplotlib artists to remove next frame
    track_lbls = []  # name annotations to remove next frame

    # Footprint plotter on near-sided
    def plot_footprint(ax, center_lon, center_lat, radius_deg):
        npts = 100
        angles = np.linspace(0, 2 * np.pi, npts)
        lons = center_lon + radius_deg * np.cos(angles)
        lats = center_lat + radius_deg * np.sin(angles)
        x, y = map2(lons, lats)
        ax.scatter(x, y, color='yellow', alpha=0.3, zorder=10)

    # ────────────────────────────────────────────────────────────────────────
    # Animate callback (called once per second)
    def animate(frame_idx, sel_dict):
        nonlocal alt_hist, frame_hist, track_objs, track_lbls

        # Clear previously drawn artists/labels
        for obj in track_objs:
            try: obj.remove()
            except Exception: pass
        track_objs.clear()
        for lbl in track_lbls:
            try: lbl.remove()
            except Exception: pass
        track_lbls.clear()

        # Reset gauges
        ax_az.cla(); init_az_compass(ax_az)
        ax_el.cla(); init_el_gauge(ax_el)

        # Reset altitude pane (ax4)
        ax4.clear()
        ax4.set_facecolor('darkslategray')
        ax4.set_title("Altitude vs. Time (rolling)", color='white')
        ax4.set_xlabel("Frame #", color='white')
        ax4.set_ylabel("Altitude (km)", color='white')
        ax4.tick_params(colors='white')

        # Time & state
        now = datetime.utcnow()
        state_dict = ConvertKepToStateVectors(sel_dict, use_skyfield=False)

        # Use the first selected sat for gauges + altitude trace
        first_name = next(iter(sel_dict))
        data = state_dict[first_name]

        # Scalars & first samples
        alt0  = data['alt_km']            # scalar (km)
        speed0 = data['speed_km_s'][0]
        lon0   = data['lons'][0]
        lat0   = data['lats'][0]

        # Compute pointing from QTH → sat
        try:
            az_deg, el_deg = az_el_from_geodetic(
                sat_lat_deg=lat0, sat_lon_deg=lon0, sat_alt_km=alt0,
                gs_lat_deg=my_lat, gs_lon_deg=my_lon, gs_h_m=0.0
            )
            print(f"[DEBUG] {first_name}: sat(lat={lat0:.3f}, lon={lon0:.3f}, alt={alt0:.1f} km)  "
                  f"QTH(lat={my_lat:.3f}, lon={my_lon:.3f})  →  Az={az_deg:.2f}°, El={el_deg:.2f}°")
        except Exception as e:
            print(f"[ERROR] Pointing calc failed: {e}")
            az_deg, el_deg = float('nan'), float('nan')

        # ---------- Compass rose (azimuth) ----------
        if not math.isnan(az_deg):
            theta = math.radians(az_deg % 360.0)
            ax_az.plot([0, theta], [0, 1.0], color='yellow', linewidth=3, zorder=5)
            ax_az.plot([theta], [1.0], marker='o', markersize=8,
                       markeredgecolor='white', markerfacecolor='yellow', zorder=6)
            ax_az.text(0.5, -0.17, f"{az_to_compass(az_deg)}  ({az_deg % 360:6.1f}°)",
                       transform=ax_az.transAxes, ha='center', va='top',
                       color='white', fontsize=11)

        # ---------- Elevation semicircle ----------
        if not math.isnan(el_deg):
            el_disp = max(0.0, min(90.0, el_deg))
            theta_el = math.radians(el_disp)    # 0..90 mapped left→up
            color_el = 'yellow' if el_deg >= 0 else 'dimgray'
            ax_el.plot([math.radians(0), theta_el], [0, 1.0], color=color_el, linewidth=3, zorder=5)
            ax_el.plot([theta_el], [1.0], marker='o', markersize=8,
                       markeredgecolor='white', markerfacecolor=color_el, zorder=6)
            # ax_el.text(0.5, -0.17, f"{el_deg:5.1f}°",
            #            transform=ax_el.transAxes, ha='center', va='top',
            #            color='white', fontsize=11)

        # ── Altitude history plot (bottom-right)
        alt_hist.append(alt0)
        frame_hist.append(frame_idx)
        ax4.plot(frame_hist, alt_hist, 'y-', linewidth=2)
        ax4.set_xlim(max(0, frame_idx - max_history + 1), frame_idx + 1)
        ax4.set_ylim(min(alt_hist) - 10, max(alt_hist) + 10)
        ax4.tick_params(colors='white')

        # ── Maps: Near-sided (ax2) and Global (ax1)
        for sat_name in sel_dict:
            sat_data = state_dict[sat_name]
            lons = sat_data['lons']
            lats = sat_data['lats']
            footprint_deg = (sat_data['alt_km'] / (6371.0 + sat_data['alt_km'])) * (180.0 / np.pi)

            # Near-sided refresh
            ax2.clear()
            map2.lon_0 = lon0
            map2.lat_0 = lat0
            map2.satellite_height = sat_data['alt_km'] * 1000.0
            map2.drawmapboundary(fill_color='aqua')
            map2.fillcontinents(color='coral', lake_color='aqua')
            map2.drawcoastlines()
            map2.drawparallels(np.arange(-90, 120, 30))
            map2.drawmeridians(np.arange(0, 420, 60))
            map2.drawstates(linewidth=1.5, color='black', zorder=10)
            map2.drawcountries(linewidth=1.5, color='black', zorder=10)

            plot_footprint(ax2, lon0, lat0, footprint_deg)

            # Ground station on near-sided
            xg, yg = map2(my_lon, my_lat)
            ax2.plot(xg, yg, 'go', markersize=8)
            ax2.annotate('Me', xy=(xg, yg), xytext=(xg + 5, yg + 5), color='white')

            # Satellite pos on near-sided
            xsn, ysn = map2(lon0, lat0)
            p2, = ax2.plot(xsn, ysn, 'r*', markersize=10, zorder=10)
            track_objs.append(p2)

            # Ground track on near-sided
            segments = split_segments(lons, lats)
            for seg_lons, seg_lats in segments:
                if len(seg_lons) < 2:
                    continue
                xx2, yy2 = map2(seg_lons, seg_lats)
                l2, = ax2.plot(xx2, yy2, color='lime', linewidth=2, zorder=5)
                d2, = ax2.plot(xx2, yy2, 'wo', markersize=2, zorder=6)
                track_objs.extend([l2, d2])

            ax2.set_xlim(map2.xmin, map2.xmax)
            ax2.set_ylim(map2.ymin, map2.ymax)
            ax2.set_title(f"Near-Sided @ {int(sat_data['alt_km'])} km", color='white')

            # GLOBAL (ax1)
            xg1, yg1 = map1(lon0, lat0)
            p1, = ax1.plot(xg1, yg1, 'r*', markersize=8, zorder=10)
            track_objs.append(p1)

            for seg_lons, seg_lats in segments:
                if len(seg_lons) < 2:
                    continue
                xx1, yy1 = map1(seg_lons, seg_lats)
                l1, = ax1.plot(xx1, yy1, color='lime', linewidth=2, zorder=5)
                d1, = ax1.plot(xx1, yy1, 'wo', markersize=2, zorder=6)
                track_objs.extend([l1, d1])

            lbl = ax1.annotate(
                sat_name,
                xy=(xg1, yg1),
                xytext=(xg1 + 6, yg1 + 6),
                color='yellow',
                fontsize=9,
                zorder=11
            )
            track_lbls.append(lbl)

        ax1.set_xlim(map1.xmin, map1.xmax)
        ax1.set_ylim(map1.ymin, map1.ymax)
        # Re-draw ground station globally
        map1.plot(x_qth, y_qth, 'go', markersize=8)
        ax1.annotate('Me', xy=(x_qth, y_qth), xytext=(x_qth + 5, y_qth + 5), color='white')

        # Title with concise status
        plt.suptitle(
            f"Amateur Satellite Tracking — UTC: {now:%Y-%m-%d %H:%M:%S}  "
            f"| {first_name}  Az {az_deg:5.1f}°  El {el_deg:5.1f}°  Alt {alt0:6.1f} km  V {speed0:4.2f} km/s",
            color='white'
        )

        fig.canvas.draw()
        fig.canvas.flush_events()
        return track_objs  # blit=False

    # Build FuncAnimation
    ani = animation.FuncAnimation(
        fig,
        animate,
        fargs=(selected,),
        interval=600,   # ms
        blit=False,
        repeat=False
    )
    plt.tight_layout()

    plt.show()

def SetupWindow(root, tle_dict, my_lat, my_lon):
    """Create a plain-white Tk window with checkboxes and two buttons."""
    root.title("Satellite Selector")
    root.configure(bg="white")

    frame1 = tk.Frame(root, bg="white")
    frame1.grid(row=0, column=0, padx=10, pady=10)

    checkbox_dict = {}
    row_ctr = 0
    col_ctr = 0

    for sat_name in sorted(tle_dict.keys()):
        var = tk.IntVar(value=0)
        cb = tk.Checkbutton(
            frame1,
            text=sat_name,
            variable=var,
            fg="black",
            bg="white",
            anchor="w"
        )
        cb.grid(row=row_ctr, column=col_ctr, sticky=tk.W, pady=2)
        checkbox_dict[sat_name] = var

        row_ctr += 1
        if row_ctr >= 20:
            row_ctr = 0
            col_ctr += 1

    frame2 = tk.Frame(root, bg="white")
    frame2.grid(row=1, column=0, padx=10, pady=10)

    run_btn = tk.Button(
        frame2,
        text="Run Prediction",
        command=lambda: runPredictionTool(checkbox_dict, tle_dict, my_lat, my_lon)
    )
    run_btn.grid(row=0, column=0, padx=5, pady=5)

    quit_btn = tk.Button(
        frame2,
        text="Quit Program",
        command=root.quit
    )
    quit_btn.grid(row=0, column=1, padx=5, pady=5)

    return checkbox_dict

if __name__ == "__main__":
    from keplerian_parser import ParseTwoLineElementFile
    tle_dict = ParseTwoLineElementFile("amateur.tle")

    # Ground station (example):
    my_lat = 41.08254088278139
    my_lon = -112.04230743099451

    # Build and show the Tk window on a white background:
    root = tk.Tk()
    checkbox_dict = SetupWindow(root, tle_dict, my_lat, my_lon)
    root.mainloop()

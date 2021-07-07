// PROCESSORS ==================================================================

function processMake(layout, selected, order, name) {
    let file = function(order, inds) {
        let f = (layout.length > 0 ? order.map(e => inds[e]) : [])
        return name(f)
    }

    switch(layout.length) {
        case 0:
            return function() {
                return { "x": 0, "y": 0, "i": [], "file": file([]) }
            }
        case 1:
            return function(A, iA) {
                return { "x": iA, "y": 0, "file": file(order, [A]) }
            }
        case 2:
            return function(A, B, iA, iB) {
                return { "x": iA, "y": iB, "file": file(order, [A, B]) }
            }
        case 3:
            return function(A, B, C, iA, iB, iC) {
                let cn = selected[layout[2]].length
                return { "x": iA*cn + iC, "y": iB, "file": file(order, [A, B, C]) }
            }
        case 4:
            return function(A, B, C, D, iA, iB, iC, iD) {
                let cn = selected[layout[2]].length
                let dn = selected[layout[3]].length
                return { "x": iA*cn + iC, "y": iB*dn + iD, "file": file(order, [A, B, C, D]) }
            }
        case 5:
            return function(A, B, C, D, E, iA, iB, iC, iD, iE) {
                let cn = selected[layout[2]].length
                let dn = selected[layout[3]].length
                let en = selected[layout[4]].length
                return { "x": iA*cn*en + iC*en + iE, "y": iB*dn + iD, "file": file(order, [A, B, C, D, E]) }
            }
        case 6:
            return function(A, B, C, D, E, F, iA, iB, iC, iD, iE, iF) {
                let cn = selected[layout[2]].length
                let dn = selected[layout[3]].length
                let en = selected[layout[4]].length
                let fn = selected[layout[5]].length
                return { "x": iA*cn*en + iC*en + iE, "y": iB*dn*fn + iD*fn + iF, "file": file(order, [A, B, C, D, E, F]) }
            }
    }
}

function processGrid(layout, selected, make) {
    switch(layout.length) {
        case 0:
            return {
                "files": files = [make()],
                "marginLeft": 5,
                "marginTop": 5,
                "nCols": 1,
                "nRows": 1
            }
        case 1:
            return {
                "files": compileSingleFiles(layout, selected, make),
                "marginLeft": 5,
                "marginTop": LABEL_SIZE,
                "nCols": selected[layout[0]].length,
                "nRows": 1
            }
        case 2:
            return {
                "files": compileDoubleFiles(layout, selected, make),
                "marginLeft": LABEL_SIZE,
                "marginTop": LABEL_SIZE,
                "nCols": selected[layout[0]].length,
                "nRows": selected[layout[1]].length
            }
        case 3:
            return {
                "files": compileTripleFiles(layout, selected, make),
                "marginLeft": LABEL_SIZE,
                "marginTop": 2*LABEL_SIZE + LABEL_PADDING,
                "nCols": selected[layout[0]].length*selected[layout[2]].length,
                "nRows": selected[layout[1]].length
            }
        case 4:
            return {
                "files": compileQuadrupleFiles(layout, selected, make),
                "marginLeft": 2*LABEL_SIZE + LABEL_PADDING,
                "marginTop": 2*LABEL_SIZE + LABEL_PADDING,
                "nCols": selected[layout[0]].length*selected[layout[2]].length,
                "nRows": selected[layout[1]].length*selected[layout[3]].length
            }
        case 5:
            return {
                "files": compileQuintupleFiles(layout, selected, make),
                "marginLeft": 2*LABEL_SIZE + LABEL_PADDING,
                "marginTop": 3*LABEL_SIZE + 2*LABEL_PADDING,
                "nCols": selected[layout[0]].length*selected[layout[2]].length*selected[layout[4]].length,
                "nRows": selected[layout[1]].length*selected[layout[3]].length
            }
        case 6:
            return {
                "files": compileSextupleFiles(layout, selected, make),
                "marginLeft": 3*LABEL_SIZE + 2*LABEL_PADDING,
                "marginTop": 3*LABEL_SIZE + 2*LABEL_PADDING,
                "nCols": selected[layout[0]].length*selected[layout[2]].length*selected[layout[4]].length,
                "nRows": selected[layout[1]].length*selected[layout[3]].length*selected[layout[5]].length
            }
    }
}

// PLOTTERS ====================================================================

function plotPath(g, S) {
    g.append("path")
        .attr("d", function(d) {
            let xscale = (d.scale ? S.xscale[d.scale.x] : S.xscale)
            let yscale = (d.scale ? S.yscale[d.scale.y] : S.yscale)
            let makePath = d3.line()
                .x(m => xscale(m))
                .y((m,i) => yscale(d.y[i]))
            return makePath(d.x)
        })
        .attr("fill", "none")
        .attr("stroke", d => (d.stroke ? d.stroke : "#555"))
        .attr("stroke-width", d => (d.width ? d.width : 1))
        .attr("opacity", d => (d.opacity ? d.opacity : null))
        .attr("stroke-linecap", d => (d.linecap ? d.linecap : null))
        .attr("stroke-dasharray", d => (d.dash ? d.dash : null))
        .attr("stroke-dashoffset", d => (d.offset ? d.offset : null))
}

function plotArea(g, S) {
    g.append("path")
        .attr("d", function(d) {
            let xscale = (d.scale ? S.xscale[d.scale.x] : S.xscale)
            let yscale = (d.scale ? S.yscale[d.scale.y] : S.yscale)
            return d3.area()
                .x(m => xscale(m))
                .y0((m, i) => yscale(d.max[i]))
                .y1((m, i) =>yscale(d.min[i]))(d.x)
        })
        .attr("fill", d => (d.fill ? d.fill : "#555"))
        .attr("stroke", d => (d.stroke ? d.stroke : "none"))
        .attr("stroke-width", d => (d.width ? d.width : "none"))
        .attr("opacity", d => (d.opacity ? d.opacity : null))
}

function plotViolin(g, S) {
    let lineForward = d3.line()
        .x(d => S.xscale(d.fx + d.dx))
        .y(d => S.yscale(d.fy + d.dy))
        .curve(d3.curveCardinal.tension(0.7))
    let lineBackward = d3.line()
        .x(d => S.xscale(d.rx + d.dx))
        .y(d => S.yscale(d.ry + d.dy))
        .curve(d3.curveCardinal.tension(0.7))

    let makeViolin = function(d) {
        if (d.direction == "horizontal") {
            var hist = d.x.map(function(e, i) { return {
                "v": e,
                "fy": e, "fx": d.y[i],
                "ry": -e, "rx": d.y[i],
                "dy": d.offset + 0.5, "dx": 0
            }})
        } else {
            var hist = d.y.map(function(e, i) { return {
                "v": e,
                "fx": e, "fy": d.x[i],
                "rx": -e, "ry": d.x[i],
                "dx": d.offset + 0.5, "dy": 0
            }})
        }

        let splits = []
        let ind = 0

        for (let i = 1; i < hist.length - 1; i++) {
            if (hist[i - 1].v == 0 && hist[i].v != 0) { ind = i }
            if (hist[i + 1].v == 0 && hist[i].v != 0) {
                let sub = hist.slice(ind - 1, i + 2)
                splits.push(lineForward(sub) +  lineBackward(sub.reverse()).replace("M", "L") + "z")
            }
        }

        return splits.join(" ")
    }

    g.append("path")
        .attr("d", d => makeViolin(d))
        .attr("fill", d => (d.fill ? d.fill : "none"))
        .attr("stroke", d => (d.stroke ? d.stroke : "#555"))
        .attr("stroke-width", d => (d.width ? d.width : 1))
        .attr("stroke-dasharray", d => (d.dash ? d.dash : null))
        .attr("opacity", d => (d.opacity ? d.opacity : null))
}

function plotSymbol(g, S) {
    g.selectAll("use")
        .data(function(d) {
            if (d.scale) {
                return d.cx.map(function(e, i) {
                    return {
                        "link": d.link[i],
                        "cx": S.xscale(e),
                        "cy": S.yscale(d.cy[i]),
                        "fill": d.fill[i],
                        "stroke": (d.stroke ? d.stroke[i] : "none"),
                        "width": (d.width ? d.width[i] : "1px"),
                    }
                })
            } else {
                if (d.bounds) { var diam = d.bounds*2 }
                else { var diam = S.axis.x.bounds[1]*2 }

                let scale = Math.min(S.subpanel.h/(diam + 1)/2, S.subpanel.w/(diam + 1)/2)
                return d.cx.map(function(e, i) {
                    return {
                        "link": d.link[i],
                        "cx": (S.subpanel.w/2 + scale*e),
                        "cy": (S.subpanel.h/2 + scale*d.cy[i]),
                        "fill": d.fill[i],
                        "stroke": d.stroke[i],
                        "width": (d.width ? d.width[i] : "1px")
                    }
                })
            }
        })
        .enter().append("use")
        .attr("transform", d => "translate(" + d.cx + "," + d.cy + ")")
        .attr("xlink:href", d => d.link)
        .attr("fill", d => d.fill)
        .attr("stroke", d => d.stroke)
        .attr("stroke-width", d => d.width)
}

function plotCircle(g, S) {
    let R = Math.min(5, Math.max(2, Math.min(S.subpanel.dw, S.subpanel.dh)/100))
    g.selectAll("circle")
        .data(function(d) {
            let xscale = (d.scale ? S.xscale[d.scale.x] : S.xscale)
            let yscale = (d.scale ? S.yscale[d.scale.y] : S.yscale)

            return d.x.map(function(e, i) {
                return {
                    "x": xscale(e),
                    "y": yscale(d.y[i]),
                    "r": (d.r ? (Array.isArray(d.r) ? d.r[i] : d.r ) : d.R ? Math.min(xscale(d.R), yscale(d.R)) : R),
                    "fill": (d.fill ? (Array.isArray(d.fill) ? d.fill[i] : d.fill) : "#555"),
                    "stroke": (d.stroke ? (Array.isArray(d.stroke) ? d.stroke[i] : d.stroke) : null),
                    "opacity": (d.opacity ? (Array.isArray(d.opacity) ? d.opacity[i] : d.opacity) : null),
                }
            })
        })
        .enter().append("circle")
            .attr("cx", d => d.x)
            .attr("cy", d => d.y)
            .attr("r", d => d.r)
            .attr("fill", d => d.fill)
            .attr("stroke", d => d.stroke)
            .attr("opacity", d => d.opacity)
}

function plotRect(g, S) {
    g.selectAll("rect")
        .data(function(d) {
            let xscale = (d.scale ? S.xscale[d.scale.x] : S.xscale)
            let yscale = (d.scale ? S.yscale[d.scale.y] : S.yscale)

            return d.x.map(function(e, i) {
                return {
                    "x": xscale(e) + (d.dx ? d.dx[i] : 0),
                    "y": yscale(d.y[i]) + (d.dy ? d.dy[i] : 0),
                    "w": xscale(d.width[i]) - xscale(0) + (d.dw ? d.dw[i] : 0),
                    "h": yscale(0) - yscale(d.height[i]) + (d.dh ? d.dh[i] : 0),
                    "fill": (d.fill ? (Array.isArray(d.fill) ? d.fill[i] : d.fill) : "#555"),
                    "stroke": (d.stroke ? (Array.isArray(d.stroke) ? d.stroke[i] : d.stroke) : null),
                    "opacity": (d.opacity ? (Array.isArray(d.opacity) ? d.opacity[i] : d.opacity) : null),
                    "dash": (d.dash ? (Array.isArray(d.dash) ? d.dash[i] : d.dash) : null),
                }
            })
        })
        .enter().append("rect")
            .attr("x", d => d.x)
            .attr("y", d => d.y)
            .attr("width", d => d.w)
            .attr("height", d => d.h)
            .attr("fill", d => d.fill)
            .attr("stroke", d => d.stroke)
            .attr("opacity", d => d.opacity)
            .attr("stroke-dasharray", d => d.dash)
}

// LABELERS ====================================================================

function labelGrid(S, P) {
    switch(S.layout.length) {
        case 0: return labelNone(S, P)
        case 1: return labelOne(S, P)
        case 2: return labelTwo(S, P)
    }
}

function labelNone(S, P) { return [] }

function labelOne(S, P) {
    let labels = []
    let layout = S.layout

    let L = layout.map(e => S.selected[e].filter(f => f != ""))

    let outerX = function(e, i) {
        return makeHorzLabel(S.panel.w, PANEL_PADDING/2 + S.panel.dw*i, 0,
            LABELS[layout[0]][e], shadeColor("#aaaaaa", i/L[0].length))
    }

    L[0].map((e, i) => labels.push(outerX(e, i)))
    return labels
}

function labelTwo(S, P) {
    let labels = []
    let layout = S.layout
    if (Array.isArray(S.layout[0])) { layout = S.selected.ordering }

    let L = layout.map(e => S.selected[e].filter(f => f != ""))

    let outerX = function(e, i) {
        return makeHorzLabel(S.panel.w, PANEL_PADDING/2 + S.panel.dw*i, 0,
            LABELS[layout[0]][e], shadeColor("#aaaaaa", i/L[0].length))
    }

    let outerY = function(e, i) {
        return makeVertLabel(S.panel.h, 0, PANEL_PADDING/2 + S.panel.dh*i,
            LABELS[layout[1]][e], shadeColor("#aaaaaa", i/L[1].length))
    }

    L[0].map((e, i) => labels.push(outerX(e, i)))
    L[1].map((e, i) => labels.push(outerY(e, i)))

    return labels
}

// DECORATORS ==================================================================

function decorateTicks(g, S, i, p) {
    addBorder(g, S.subpanel.w, S.subpanel.h, "#ccc")

    // Create and align groups for ticks.
    let dx = alignHorzAxis(S, i)
    let dy = alignVertAxis(S, i)

    // Create group to hold ticks.
    let G = S.G.append("g")
        .attr("id", "ticks")
        .attr("transform", "translate(" + dx + "," + dy + ")")

    let A = S.axis

    // Create ticks.
    let ticks = []
    ticks.push(makeHorzTicks(S, 0, S.subpanel.h, A.x))
    ticks.push(makeVertTicks(S, 0, 0, A.y))

    // Create axis labels.
    let labels = []
    labels.push(makeHorzLabel(S.subpanel.w, 0, alignHorzText(S), A.x.title, "none"))
    labels.push(makeVertLabel(S.subpanel.h, alignVertText(S), 0, A.y.title, "none"))

    addTicks(G, ticks)
    addLabels(G, labels)
}

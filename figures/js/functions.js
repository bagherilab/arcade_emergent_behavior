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

function plotSymbol(g, S) {
    g.selectAll("use")
        .data(function(d) {
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
        })
        .enter().append("use")
        .attr("transform", d => "translate(" + d.cx + "," + d.cy + ")")
        .attr("xlink:href", d => d.link)
        .attr("fill", d => d.fill)
        .attr("stroke", d => d.stroke)
        .attr("stroke-width", d => d.width)
}

// LABELERS ====================================================================

function labelGrid(S, P) {
    switch(S.layout.length) {
        case 0: return labelNone(S, P)
        case 2: return labelTwo(S, P)
    }
}

function labelNone(S, P) { return [] }

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


function saveAs(filename) {
    var filename = filename.split("/")
    filename = filename[filename.length - 1]
    let svg = document.getElementById("download")
    let serializer = new XMLSerializer()
    let source = serializer.serializeToString(svg)
    let svgBlob = new Blob([source], {type:"image/svg+xml;charset=utf-8"})
    let downloadLink = document.createElement("a")
    downloadLink.href = URL.createObjectURL(svgBlob)
    downloadLink.download = filename + ".svg"
    document.body.appendChild(downloadLink)
    downloadLink.click()
    document.body.removeChild(downloadLink)
}

// -----------------------------------------------------------------------------

function shadeColor(color, percent) {
    if (color.slice(0,3) == "rgb") {
        let split = color.split(",")
        var R = Number(split[0].split("(")[1])
        var G = Number(split[1])
        var B = Number(split[2].split(")")[0])
    }
    else {
        let f = parseInt(color.slice(1), 16)
        var R = f >> 16
        var G = f >> 8&0x00FF
        var B = f&0x0000FF
    }

    let t = percent < 0 ? 0:255
    let p = percent < 0 ? percent*-1:percent

    let r = (Math.round((t - R)*p)+R)*0x10000
    let g = (Math.round((t - G)*p)+G)*0x100
    let b = (Math.round((t - B)*p)+B)

    return "#" + (0x1000000 + r + g + b).toString(16).slice(1)
}

function linspace(start, end, n) {
    let delta = (end - start)/(n - 1)
    let bins = d3.range(start, end + delta, delta).slice(0, n)
    return bins.map(e => Number((parseFloat(e).toPrecision(12))))
}

function linrange(start, n) {
    let r = []
    for (let i = 0; i < n; i++) { r.push(start + i) }
    return r
}

function contains(val, arr) {
    for (let i = 0; i < arr.length; i++) { if (arr[i] == val) { return true }}
    return false
}

// -----------------------------------------------------------------------------

function findNaNs(arr) {
    return arr.map((e, i) => !(Number.isNaN(e) || e == "nan"))
}

function removeNaNs(arr, remove) {
    return arr.filter((e, i) => remove[i])
}

function binData(data, bounds, bandwidth) {
    let N = data.length*bandwidth
    let n = (bounds[1] - bounds[0])/bandwidth
    let lower = bounds[0] - bandwidth/2
    let upper = bounds[1] + bandwidth/2
    let bins = linspace(lower, upper - bandwidth, n + 1)

    let hist = d3.histogram().domain([lower, upper]).thresholds(bins)(data)
        .map(e => ({
            "n": e.length,
            "x": e.x0,
            "y": e.length/N,
            "m": Number(((e.x0 + e.x1)/2).toFixed(3)),
        }))

    return hist
}

// -----------------------------------------------------------------------------

function makeHex() {
    let points = [0, 1, 2, 3, 4, 5, 0]
    let theta = points.map(e => Math.PI*(60*e)/180.0)
    let dx = theta.map(e => 2/Math.sqrt(3)*Math.cos(e))
    let dy = theta.map(e => 2/Math.sqrt(3)*Math.sin(e))
    return points.map((e, i) => [dx[i], dy[i]])
}

function makeTri(i, n) {
    let hex = makeHex()
    let tri = []
    for (let p = i; p <= i + n; p++) { tri.push([hex[p%6][0], hex[p%6][1]]) }
    return [[0,0]].concat(tri)
}

// -----------------------------------------------------------------------------

function makeVertLabel(h, x, y, text, fill) {
    return {
        "w": LABEL_SIZE, "h": h, "x": x - LABEL_SIZE, "y": y, "text": text, "fill": fill,
        "tx": x - LABEL_SIZE + FONT_SIZE + FONT_PADDING, "ty": y + h/2, "rotate": true
    }
}

function makeHorzLabel(w, x, y, text, fill) {
    return {
        "w": w, "h": LABEL_SIZE, "x": x, "y": y - LABEL_SIZE, "text": text, "fill": fill,
        "tx": x + w/2, "ty": y - FONT_PADDING, "rotate": false
    }
}

// -----------------------------------------------------------------------------

function makeVertTicks(S, x, y, axis, scale) {
    let bounds = axis.bounds
    let padding = axis.padding ? axis.padding : 0
    let ticks = linspace(bounds[0], bounds[1], axis.n).map(d => (Math.abs(d) < 10E-10 ? 0 : d))
    let t = []

    if (!scale) {
        scale = d3.scaleLinear().range([S.panel.h - SUBPANEL_PADDING - AXIS_PADDING.bottom, 0]).domain([bounds[0] - padding, bounds[1] + padding])
    }

    for (let i = 0; i < ticks.length; i++) {
        t.push({
            "tx": x - (FONT_SIZE - 2) - 3 - FONT_PADDING,
            "ty": y + scale(ticks[i]) + (FONT_SIZE - 2)/2,
            "y1": y + scale(ticks[i]),
            "y2": y + scale(ticks[i]),
            "x1": x,
            "x2": x - 3,
            "text": (axis.labels ? axis.labels(ticks[i], i) : ticks[i]),
        })
    }

    return t
}

function makeHorzTicks(S, x, y, axis, scale) {
    let bounds = axis.bounds
    let padding = axis.padding ? axis.padding : 0
    let ticks = linspace(bounds[0], bounds[1], axis.n).map(d =>(Math.abs(d) < 10E-10 ? 0 : d))
    let t = []

    if (!scale) {
        scale = d3.scaleLinear().range([0, S.subpanel.w]).domain([bounds[0] - padding, bounds[1] + padding])
    }

    for (let i = 0; i < ticks.length; i++) {
        t.push({
            "tx": x + scale(ticks[i]),
            "ty": y + 3 + FONT_SIZE - 2 + FONT_PADDING,
            "y1": y,
            "y2": y + 3,
            "x1": x + scale(ticks[i]),
            "x2": x + scale(ticks[i]),
            "text": (axis.labels ? axis.labels(ticks[i], i) : ticks[i]),
        })
    }

    return t
}

// -----------------------------------------------------------------------------

function alignVertAxis(S, i) {
    return PANEL_PADDING/2 + SUBPANEL_PADDING/2 + S.margin.axis.top
        + (i.length < 2 ? 0 : i.length < 4 ? S.panel.dh*i[1] :
            S.panel.dh*i[1]*S.selected[S.layout[3]].length + S.panel.dh*i[3])
}

function alignHorzAxis(S, i) {
    return PANEL_PADDING/2 + SUBPANEL_PADDING/2 + S.margin.axis.left
        + (i.length == 0 ? 0 : i.length < 3 ? S.panel.dw*i[0] :
            S.panel.dw*i[0]*S.selected[S.layout[2]].length + S.panel.dw*i[2])
}

function alignVertText(S) {
    return -5 - 2*FONT_PADDING - (FONT_SIZE - 2)*2
}

function alignHorzText(S) {
    return S.subpanel.h + LABEL_SIZE + 5 + FONT_SIZE - 2 + 2*FONT_PADDING
}

// -----------------------------------------------------------------------------

function compileSingleFiles(layout, selected, make) {
    let files = []
    let len = layout.map(d => selected[d].length)
    let A = selected[layout[0]]

    for (let a = 0; a < len[0]; a++) {
        let f = make(A[a], a)
        files.push({
            "file": f.file,
            "x": f.x,
            "y": f.y,
            "i": [a]
        })
    }

    return files
}

function compileDoubleFiles(layout, selected, make) {
    let files = []
    let len = layout.map(d => selected[d].length)
    let A = selected[layout[0]]
    let B = selected[layout[1]]

    for (let a = 0; a < len[0]; a++) {
        for (let b = 0; b < len[1]; b++) {
            let ab = make(A[a], B[b], a, b)
            files.push({
                "file": ab.file,
                "x": ab.x,
                "y": ab.y,
                "i": [a, b]
            })
        }
    }

    return files
}

function compileTripleFiles(layout, selected, make) {
    let files = []
    let len = layout.map(d => selected[d].length)
    let A = selected[layout[0]]
    let B = selected[layout[1]]
    let C = selected[layout[2]]

    for (let a = 0; a < len[0]; a++) {
        for (let b = 0; b < len[1]; b++) {
            for (let c = 0; c < len[2]; c++) {
                let abc = make(A[a], B[b], C[c], a, b, c)
                files.push({
                    "file": abc.file,
                    "x": abc.x,
                    "y": abc.y,
                    "i": [a, b, c] })
            }
        }
    }

    return files
}

// -----------------------------------------------------------------------------

function addBorder(g, width, height, stroke) {
    return g.append("rect")
        .attr("width", width)
        .attr("height", height)
        .attr("fill", "none")
        .attr("stroke", stroke)
        .attr("stroke-width", "1px")
}

function addLine(g, x1, x2, y1, y2, stroke, width) {
    return g.append("line")
        .attr("x1", x1)
        .attr("x2", x2)
        .attr("y1", y1)
        .attr("y2", y2)
        .attr("stroke", stroke)
        .attr("stroke-width", width)
}

function addLabels(g, labels) {
    let G = g.append("g").attr("id", "labels")

    G.selectAll("rect").data(labels)
        .enter().append("rect")
        .attr("width", d => d.w)
        .attr("height", d => d.h)
        .attr("x", d => d.x)
        .attr("y", d => d.y)
        .attr("fill", d => d.fill)

    G.selectAll("text").data(labels)
        .enter().append("text")
        .html(d => d.text)
        .attr("transform", d => d.rotate ? "rotate(-90," + d.tx + "," + d.ty + ")" : null)
        .attr("font-size", FONT_SIZE + "pt")
        .attr("font-weight", "bold")
        .attr("font-family", "Helvetica")
        .attr("text-anchor", "middle")
        .attr("x", d => d.tx)
        .attr("y", d => d.ty)
}

function addTicks(g, ticks) {
    let G = g.selectAll("g").data(ticks).enter().append("g")

    G.selectAll("line")
        .data(d => d).enter().append("line")
            .attr("x1", d => d.x1 )
            .attr("x2", d => d.x2 )
            .attr("y1", d => d.y1 )
            .attr("y2", d => d.y2 )
            .attr("stroke", "#000")

    G.selectAll("text")
        .data(d => d).enter().append("text")
            .html(d => d.text)
            .attr("font-size", (FONT_SIZE - 2) + "pt")
            .attr("font-family", "Helvetica")
            .attr("text-anchor", "middle")
            .attr("x", d => d.tx)
            .attr("y", d => d.ty)
}

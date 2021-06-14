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

// -----------------------------------------------------------------------------

function addBorder(g, width, height, stroke) {
    return g.append("rect")
        .attr("width", width)
        .attr("height", height)
        .attr("fill", "none")
        .attr("stroke", stroke)
        .attr("stroke-width", "1px")
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


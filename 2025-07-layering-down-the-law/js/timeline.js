// Timeline data
const editions = [
    { year: 1985, label: "1st Edition" },
    { year: 1988, label: "2nd Edition" },
    { year: 1992, label: "3rd Edition" },
    { year: 1996, label: "4th Edition" },
    { year: 2001, label: "5th Edition" },
    // { year: 2024, label: "12th Edition" },
    { year: 2003, label: "Amy reads the 5th edition", above: true }
];

const notes = [
];

// Timeline configuration
const margin = { top: 40, right: 50, bottom: 60, left: 50 };
const width = 1600 - margin.left - margin.right;
const height = 900 - margin.bottom - margin.top;

// Create timeline visualization
function createTimeline() {
    // Create SVG
    const svg = d3.select("#timeline-chart")
        .append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.bottom + margin.top)
        .style("margin-top", "2em");

    const g = svg.append("g")
        .attr("transform", `translate(${margin.left},${margin.top})`);

    // Create scale
    const xScale = d3.scaleLinear()
        // .domain([1693, 2025])
        .domain([1980, 2025])
        .range([0, width]);

    // Create and call axis
    const xAxis = d3.axisTop(xScale)
        .ticks(20)
        .tickFormat(d3.format("d"))
        .tickSize(15);

    const axisGroup = g.append("g")
        .attr("transform", `translate(0,${height/2})`)
        .call(xAxis);
    
    axisGroup.selectAll("path")
        .attr("stroke-width", "1px");
    
    axisGroup.selectAll("text")
        .attr("font-size", "10px");

        // Create and call axis
    const xAxisMinor = d3.axisTop(xScale)
        .ticks(10)
        .tickFormat("")
        .tickSize(6);

    const axisGroupMinor = g.append("g")
        .attr("transform", `translate(0,${height/2})`)
        .call(xAxisMinor);
    
    axisGroupMinor.selectAll("path")
        .attr("stroke-width", "1px");
    
    axisGroupMinor.selectAll("text")
        .attr("font-size", "10px");

    // Edition markers
    const editionGroup = g.append("g")
        .attr("class", "editions");

    editionGroup.selectAll("circle")
        .data(editions)
        .enter()
        .append("circle")
        .attr("cx", d => xScale(d.year))
        .attr("cy", height/2)
        .attr("r", 3);

    editionGroup.selectAll("text")
        .data(editions)
        .enter()
        .append("text")
        .attr("x", d => xScale(d.year))
        .attr("y", d => d.above ? height/2 - 35 : height/2 + 15)
        .attr("text-anchor", "middle")
        .attr("font-size", "8px")
        .attr("font-weight", "bold")
        .text(d => `${d.label} (${d.year})`);
}

// Initialize timeline when DOM is ready
document.addEventListener('DOMContentLoaded', function() {
    createTimeline();
}); 
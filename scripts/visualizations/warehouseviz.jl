
using Luxor, Colors

tempHueList = collect(Colors.color_names)

#Exclude colors that are too light
hueList = []
for item in tempHueList
	if (item[2][1] + item[2][2] + item[2][3] <= 515) & (max(abs(item[2][1]-item[2][2]), abs(item[2][1]-item[2][3]), abs(item[2][2]-item[2][3])) >= 100)
		push!(hueList, item) 
	end
end

#---------------------------------------------------------------------------------------#

function warehouseviz(wsdrawingname, vizx)

	vizy = convert(Int, round(vizx * warehouse_y_length_meters / warehouse_x_length_meters, digits=0))

	#Find coordinates for each workstation and pod storage location
	locPoints = Dict()
	for l in union(storagelocs, workstations)
		newx = loccoords[l,1] * (vizx - 100) / warehouse_x_length_meters - (vizx - 100) / 2
		newy = loccoords[l,2] * (vizy - 100) / warehouse_y_length_meters - (vizy - 100) / 2
		locPoints[l] = Point((newx, newy))
	end

	#-------------------------------------------------------------------------#
	
	#Create a list of boxes for locations
	locationsquares = []
	for l in storagelocs
		corner = locPoints[l]
		width = 6 * (vizx - 100) / warehouse_x_length_meters
		height = 4 * (vizy - 100) / warehouse_y_length_meters
		thickness = 5
		push!(locationsquares, (corner, width, height, thickness))
	end
	for l in workstations
		corner = locPoints[l]
		width = 6 * (vizx - 100) / warehouse_x_length_meters
		height = 4 * (vizy - 100) / warehouse_y_length_meters
		thickness = 8
		push!(locationsquares, (corner, width, height, thickness))
	end

	intersectionmarks, intersectionlabels = [], []
	thickness = 3
	intersectioncolor = (150,0,0)
	#maxviol = maximum(values(maxtraffic))
	for l in intersections
		newx = intcoords[l][1] * (vizx - 100) / warehouse_x_length_meters - (vizx - 100) / 2
		newy = intcoords[l][2] * (vizy - 100) / warehouse_y_length_meters - (vizy - 100) / 2
		
		#Standard viz
		#line1_s = Point((newx, newy)) + Point(5*vizx/1000,5*vizx/1000)
		#line1_e = Point((newx, newy)) + Point(-5*vizx/1000,-5*vizx/1000)
		#push!(intersectionmarks, (line1_s, line1_e, thickness, intersectioncolor))
		#line2_s = Point((newx, newy)) + Point(-5*vizx/1000,5*vizx/1000)
		#line2_e = Point((newx, newy)) + Point(5*vizx/1000,-5*vizx/1000)
		#push!(intersectionmarks, (line2_s, line2_e, thickness, intersectioncolor))

		#Congestion viz
		trafficcolor = (0,0,0)
		#if intersectionviolations[l] > 1e-4 
		#	trafficcolor = ((maxtraffic[l])/(maxviol)*255,0,0)
		#else
		#	trafficcolor = ((1-maxtraffic[l])*255,(1-maxtraffic[l])*255,(1-maxtraffic[l])*255)
		#end
		push!(intersectionmarks, (Point((newx, newy)), 15*vizx/1000, thickness, trafficcolor))
		
		#Label
		push!(intersectionlabels, (Point((newx, newy)), string(l)))

	end

	#-------------------------------------------------------------------------#

	Drawing(vizx, vizy, wsdrawingname)
	origin()
	background("white")

	#Warehouse outline
	setdash("solid")
	setcolor(convert(Colors.HSV, Colors.RGB(1,1,1)))
	Luxor.rect( Point(- (vizx - 100) / 2, - (vizy - 100) / 2), vizx - 100, vizy - 100, :fill)	
	setcolor("black")
	setline(5)
	Luxor.rect( Point(- (vizx - 100) / 2, - (vizy - 100) / 2), vizx - 100, vizy - 100, :stroke)	

	#Location boxes
	for box in locationsquares
		setline(box[4])
		setcolor("white")
		Luxor.rect(box[1], box[2], box[3], :fill)	
		setcolor("black")
		Luxor.rect(box[1], box[2], box[3], :stroke)	
	end

	#Location boxes
	#for segment in intersectionmarks
	#	setline(segment[3]*vizx/1000)
	#	r_val, g_val, b_val = segment[4]
	#	setcolor(convert(Colors.HSV, Colors.RGB(r_val / 255, g_val / 255, b_val / 255)))
	#	Luxor.line(segment[1], segment[2], :stroke)
	#end

	setcolor(convert(Colors.HSV, Colors.RGB(150 / 255, 150 / 255, 150 / 255)))
	setline(30)
	toppoint, bottompoint = (horstreets[1,1]+horstreets[1,2])/2, (horstreets[size(horstreets)[1],1]+horstreets[size(horstreets)[1],2])/2
	queuepoint1, queuepoint2 = maximum(loccoords[1:workstations[1]-1,2])+4, maximum(loccoords[:,2])
	queuepoint = (queuepoint1+queuepoint2)/2
	for st in 1:size(vertstreets)[1]
		newx = (vertstreets[st,1]+vertstreets[st,2])/2 * (vizx - 100) / warehouse_x_length_meters - (vizx - 100) / 2
		y1 = toppoint * (vizy - 100) / warehouse_y_length_meters - (vizy - 100) / 2
		y2 = bottompoint * (vizy - 100) / warehouse_y_length_meters - (vizy - 100) / 2
		Luxor.line(Point(newx,y1),Point(newx,y2), :stroke)
	end
	leftpoint, rightpoint = (vertstreets[1,1]+vertstreets[1,2])/2, (vertstreets[size(vertstreets)[1],1]+vertstreets[size(vertstreets)[1],2])/2
	for st in 1:size(horstreets)[1]
		newy = (horstreets[st,1]+horstreets[st,2])/2 * (vizy - 100) / warehouse_y_length_meters - (vizy - 100) / 2
		x1 = leftpoint * (vizx - 100) / warehouse_x_length_meters - (vizx - 100) / 2
		x2 = rightpoint * (vizx - 100) / warehouse_x_length_meters - (vizx - 100) / 2
		Luxor.line(Point(x1,newy),Point(x2,newy), :stroke)
	end
	for q in queueintersections
		newx = intcoords[q][1] * (vizx - 100) / warehouse_x_length_meters - (vizx - 100) / 2
		y1 = queuepoint * (vizy - 100) / warehouse_y_length_meters - (vizy - 100) / 2
		y2 = intcoords[q][2] * (vizy - 100) / warehouse_y_length_meters - (vizy - 100) / 2
		Luxor.line(Point(newx,y1),Point(newx,y2), :stroke)
	end
	
	for segment in intersectionmarks
		setline(10)
		r_val, g_val, b_val = segment[4]
		setcolor(convert(Colors.HSV, Colors.RGB(r_val / 255, g_val / 255, b_val / 255)))
		Luxor.circle(segment[1], segment[2], :fill)

		setline(5)
		setcolor("black")
		Luxor.circle(segment[1], segment[2], :stroke)
	end

	setcolor("black")
	fontsize(100*vizx/2000)
	#for w in union(storagelocs,workstations)
	#	Luxor.text(string(w), locPoints[w] + Point(3 * (vizx - 100) / warehouse_x_length_meters, 2 * (vizy - 100) / warehouse_y_length_meters), halign=:center,   valign = :middle)
	#end
	for w in workstations[2]
		Luxor.text("w", locPoints[w] + Point(3 * (vizx - 100) / warehouse_x_length_meters, 2 * (vizy - 100) / warehouse_y_length_meters), halign=:center,   valign = :middle)
	end
	#fontsize(20*vizx/2000)
	#for lbl in intersectionlabels
	#	Luxor.text(lbl[2], lbl[1])
	#end

	finish()
	preview()

end
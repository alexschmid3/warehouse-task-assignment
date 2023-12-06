
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
	for l in intersections
		newx = intcoords[l][1] * (vizx - 100) / warehouse_x_length_meters - (vizx - 100) / 2
		newy = intcoords[l][2] * (vizy - 100) / warehouse_y_length_meters - (vizy - 100) / 2
		
		line1_s = Point((newx, newy)) + Point(5,5)
		line1_e = Point((newx, newy)) + Point(-5,-5)
		push!(intersectionmarks, (line1_s, line1_e, thickness, intersectioncolor))

		line2_s = Point((newx, newy)) + Point(-5,5)
		line2_e = Point((newx, newy)) + Point(5,-5)
		push!(intersectionmarks, (line2_s, line2_e, thickness, intersectioncolor))

		push!(intersectionlabels, (Point((newx, newy)), string(l)))

	end

	#-------------------------------------------------------------------------#

	Drawing(vizx, vizy, wsdrawingname)
	origin()
	background("white")

	#Warehouse outline
	setcolor("black")
	setdash("solid")
	setline(2)
	Luxor.rect( Point(- (vizx - 100) / 2, - (vizy - 100) / 2), vizx - 100, vizy - 100, :stroke)	

	#Location boxes
	for box in locationsquares
		setline(box[4])
		Luxor.rect(box[1], box[2], box[3], :stroke)	
	end

	#Location boxes
	for segment in intersectionmarks
		setline(segment[3])
		r_val, g_val, b_val = segment[4]
		setcolor(convert(Colors.HSV, Colors.RGB(r_val / 255, g_val / 255, b_val / 255)))
		Luxor.line(segment[1], segment[2], :stroke)
		
	end

	setcolor("black")
	fontsize(36)
	for w in union(storagelocs,workstations)
		Luxor.text(string(w), locPoints[w] + Point(3 * (vizx - 100) / warehouse_x_length_meters, 2 * (vizy - 100) / warehouse_y_length_meters), halign=:center,   valign = :middle)
	end
	#fontsize(20)
	#for lbl in intersectionlabels
	#	label(lbl[2], :N , lbl[1] + Point(0,-3))
	#end

	finish()
	preview()

end
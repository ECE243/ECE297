/* 
 * Copyright 2020 University of Toronto
 *
 * Permission is hereby granted, to use this software and associated 
 * documentation files (the "Software") in course work at the University 
 * of Toronto, or for personal use. Other uses are prohibited, in 
 * particular the distribution of the Software either publicly or to third 
 * parties.
 *
 * The above copyright notice and this permission notice shall be included in 
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#include "m1.h"
#include "StreetsDatabaseAPI.h"
#include "OSMDatabaseAPI.h"
#include <vector>
#include <unordered_set>
#include <map>

//*********************************************************************
//***********        Global Variables Declaration           ***********
//*********************************************************************

// 2D array of all street segments IDs on the map, indexed by street IDs
std::vector<std::vector<StreetSegmentIndex>> streetSegmentsOfStreets;

// 2D array of all intersections IDs on the map, indexed by street IDs
//
// The first contains a vector as its inner array (useful when accessing
// all intersections of a street at once)
std::vector<std::vector<IntersectionIndex>> intersectionsAlongStreets_v;
//
// The second contains an unordered_set as its inner array (useful for searching
// for a single intersection from the list of all intersections of a street)
std::vector<std::unordered_set<IntersectionIndex>> intersectionsAlongStreets_uos;

// Unordered map of all node_ids (the OSMDatabaseAPI indices), indexed by the 
// corresponding OSMIDs of those OSMNode objects
std::unordered_map<OSMID, int> nodeIDToOSMIDMap;
// Unordered map of all way_ids (the OSMDatabaseAPI indices), indexed by the 
// corresponding OSMIDs of those OSMWay objects
std::unordered_map<OSMID, int> wayIDToOSMIDMap;

// Vector of vectors of ints of all street segments attached to an intersection
// vector[intersectionid] = vector of streetSegments
std::vector<std::vector<StreetSegmentIndex>> intersectionStreetSegments;

// Vector storing the travel time for each street segment (indexed by the 
// street segment id)
std::vector<double> total_travel_time;

// Vector of vectors of strings of all street segments names attached to an intersection
// vector[intersectionid] = vector of streetSegmentNames
std::vector<std::vector<std::string>> intersectionStreetNames;

// Vector of vectors of ints of all adjacent intersection ids connected to an intersection
// that can be traveled down in the correct direction
// vector[intersectionid] = vector of adjacent intersections
std::vector<std::vector<IntersectionIndex>> adjacentIntersections;

// Multimap storing pairs of name and corresponding intersection ID
std::multimap<std::string, int> streetNamesDatabase;

/*prototyping the functions*/
void load_street_segments_of_streets(int numStreets);
void load_intersections_along_streets(int numStreets);
void load_node_id();
void load_way_id();
void load_street_segments_and_street_names_of_intersection(int numIntersections);
double street_length(int street_segment_id );
void load_time_taken(int numStreetSegments);
void load_adjacent_intersections(int numIntersections);
std::pair<double, double> getEquirectangularProjection(LatLon point, double avgLatitude);
double findAverageLatitude(std::vector<LatLon> points);












// ******** MOVE THIS LATER *******************************
void removeSpacesAndCaptialsInString(std::string& str);
// ********************************************************

//*********************************************************************
//***********          Load Map Helper Functions            ***********
//*********************************************************************

// Load the streetSegmentsOfStreets data structure
void load_street_segments_of_streets(int numStreets){
    // The outer array is of size 'numStreets'
    streetSegmentsOfStreets = std::vector<std::vector<StreetSegmentIndex>>(numStreets);
    
    // The i'th element of the outer array is an array of all street segment 
    // IDs that are part of the i'th street.
    //
    // To populate these inner arrays, iterate through every street segment on 
    // the map and add its ID to the streetSegmentsOfStreets 2D array under the 
    // appropriate street
    int numStreetSegments = getNumStreetSegments();
    
    for (StreetSegmentIndex i_streetSeg = 0; i_streetSeg < numStreetSegments; i_streetSeg++) {
        StreetIndex i_street = getInfoStreetSegment(i_streetSeg).streetID;
        streetSegmentsOfStreets[i_street].push_back(i_streetSeg);
    }
}

// Load the intersectionsAlongStreets data structures 
// (both the _v and _uos versions)
void load_intersections_along_streets(int numStreets){
    // The outer array is of size 'numStreets'
    intersectionsAlongStreets_v = std::vector<std::vector<IntersectionIndex>>(numStreets);
    intersectionsAlongStreets_uos = std::vector<std::unordered_set<IntersectionIndex>>(numStreets);
    
    // The i'th element of the outer array is an array of all intersection IDs 
    // that are part of the i'th street.
    //
    // To populate these inner arrays, iterate through every street on 
    // the map. 
    for (StreetIndex i_street = 0; i_street < numStreets; i_street++) {        
        // For each street, look at all the street segments that are part of it, 
        // and use them to determine all intersections that the street passes through
        for (int nSeg = 0; nSeg < streetSegmentsOfStreets[i_street].size(); nSeg++) {
            // Get the street segment object of the n'th street segment of street i
            StreetSegmentIndex i_nthStreetSeg = streetSegmentsOfStreets[i_street][nSeg];
            InfoStreetSegment i_streetSegInfo = getInfoStreetSegment(i_nthStreetSeg);
            
            // Add both 'to' and 'from' intersections of the n'th street segment 
            // to the list of intersections of street i.
            //
            // Note: to avoid duplicates, need to check if the 'to' and 'from' 
            // have already been included in the list. To achieve this quickly,
            // use the unordered_set.
            if (intersectionsAlongStreets_uos[i_street].find(i_streetSegInfo.from) == intersectionsAlongStreets_uos[i_street].end()) {
                intersectionsAlongStreets_uos[i_street].insert(i_streetSegInfo.from);
                intersectionsAlongStreets_v[i_street].push_back(i_streetSegInfo.from);
            }
            if (intersectionsAlongStreets_uos[i_street].find(i_streetSegInfo.to) == intersectionsAlongStreets_uos[i_street].end()) {
                intersectionsAlongStreets_uos[i_street].insert(i_streetSegInfo.to);
                intersectionsAlongStreets_v[i_street].push_back(i_streetSegInfo.to);
            }
        }
    }
}

// Load the nodeIDToOSMIDMap data structure
void load_node_id(){
    // Iterate thorough all nodes on the map, find their nodeIDs, and add the
    // pair (OSMID, node_id) to the nodeIDToOSMIDMap
    int numNodes = getNumberOfNodes();    
    for (int i_node = 0; i_node < numNodes; i_node++) {
        // Get the OSMID of node i
        OSMID nodeOSMID = getNodeByIndex(i_node)->id();
        
        // Add a new entry to nodeIDToOSMIDMap with
        // key = nodeOSMID and value = i_node
        nodeIDToOSMIDMap[nodeOSMID] = i_node;
    }
}

// Load the wayIDToOSMIDMap data structure
void load_way_id(){
    // Iterate thorough all ways on the map, find their wayIDs, and add the
    // pair (OSMID, way_id) to the wayIDToOSMIDMap
    int numWays = getNumberOfWays();    
    for (int i_way = 0; i_way < numWays; i_way++) {
        // Get the OSMID of way i
        OSMID wayOSMID = getWayByIndex(i_way)->id();
        
        // Add a new entry to wayIDToOSMIDMap with
        // key = wayOSMID and value = i_way
        wayIDToOSMIDMap[wayOSMID] = i_way;
    }
}

// Load the intersectionStreetSegments and Load the intersectionStreetNames
void load_street_segments_and_street_names_of_intersection(int numIntersections){
    StreetSegmentIndex currentRoad;
    intersectionStreetSegments.resize(getNumIntersections());

    InfoStreetSegment currentRoadInfo;
    std::string currentRoadName;
    intersectionStreetNames.resize(getNumIntersections());

    for(IntersectionIndex intersectionID=0; intersectionID < numIntersections; intersectionID++){
        for(int segmentNumber=0; segmentNumber < getIntersectionStreetSegmentCount(intersectionID); segmentNumber++){
            
            currentRoad = getIntersectionStreetSegment(intersectionID, segmentNumber);
            intersectionStreetSegments[intersectionID].push_back(currentRoad);
            
            currentRoadInfo = getInfoStreetSegment(currentRoad);
            currentRoadName = getStreetName(currentRoadInfo.streetID);
            intersectionStreetNames[intersectionID].push_back(currentRoadName);
        }
    }
}
double street_length(int street_segment_id ){// returns the street length after conversion to km from m
    
    double length_of_street=find_street_segment_length(street_segment_id)*3.6;
    return length_of_street;
}

void load_time_taken(int numStreetSegments){// uses speed = distance/time to find time

    InfoStreetSegment onWhichRoad; 
    total_travel_time.resize(numStreetSegments);
    for(StreetSegmentIndex streetSegment= 0 ; streetSegment < numStreetSegments; streetSegment++  ){
        onWhichRoad = getInfoStreetSegment(streetSegment); //;
        total_travel_time[streetSegment]=street_length(streetSegment)/onWhichRoad.speedLimit;

    }

}

// Load the adjacentIntersections
void load_adjacent_intersections(int numIntersections){
    std::vector<int> intersectionStreetSegments1;
    std::pair<int, int> intersectionPair;
    adjacentIntersections.resize(getNumIntersections());
    std::unordered_set<int> adjacentIntersectionsCheck; //create an unordered set to search for duplicates of order 1

    for(IntersectionIndex intersectionID=0; intersectionID < numIntersections; intersectionID++){
        intersectionStreetSegments1 = find_street_segments_of_intersection(intersectionID);
        intersectionPair.first = intersectionID;
        adjacentIntersectionsCheck.clear(); 
        
        for(int index=0; index < intersectionStreetSegments1.size(); index++){
            
            //create a vector of 2 intersection IDs, but make sure they do not repeat an intersection
            if(getInfoStreetSegment(intersectionStreetSegments1[index]).to == intersectionID){
                intersectionPair.second = getInfoStreetSegment(intersectionStreetSegments1[index]).from;
            } else {
                intersectionPair.second = getInfoStreetSegment(intersectionStreetSegments1[index]).to;
            }
            if(are_directly_connected(intersectionPair)){
                //compare the returned iterator to end of vector to see if duplicate
                if (adjacentIntersectionsCheck.find(intersectionPair.second) == adjacentIntersectionsCheck.end()) { 
                    adjacentIntersections[intersectionID].push_back(intersectionPair.second);
                    adjacentIntersectionsCheck.insert(intersectionPair.second); 
                }
            }
        }
    }    
}

#define STREETS_BIN_PATH_EXTENSION_SIZE 11 // Number of characters in "streets.bin"

bool load_map(std::string map_path) { 
    // Check if the map_path is valid based on the size (if it isn't at least
    // STREETS_BIN_PATH_EXTENSION_SIZE characters long, it is definitely invalid
    // since it means there isn't even enough space for the "streets.bin" extension)
    if (map_path.size() < STREETS_BIN_PATH_EXTENSION_SIZE) {
        return false;
    }
    
    // Using the map_path, determine the OSM path (simply change the extension
    // by removing the "streets.bin" from the end and adding "osm.bin")
    std::string osm_bin_path = map_path;
    // Remove "streets.bin" (we know there are enough elements to remove)
    osm_bin_path.resize(map_path.size() - STREETS_BIN_PATH_EXTENSION_SIZE);
    osm_bin_path.append("osm.bin"); // Append "osm.bin"
    
    // Initialize the StreetsDatabaseAPI and the OSMDatabaseAPI.
    // Ensure that they have been loaded successfully. If there is a failure,
    // exit the loop
    if (!loadStreetsDatabaseBIN(map_path) || !loadOSMDatabaseBIN(osm_bin_path)) {
        return false;
    }
    
    // Loading all map related data structures
    int numStreets = getNumStreets();
    int numIntersections = getNumIntersections();
    int numStreetSegments = getNumStreetSegments();
    
    // Load the streetSegmentsOfStreets data structure
    load_street_segments_of_streets(numStreets);
    
    // Load the intersectionsAlongStreets data structures 
    // (both the _v and _uos versions)
    load_intersections_along_streets(numStreets);
    
    // Load the nodeIDToOSMIDMap data structure
    load_node_id();
    
    // Load the wayIDToOSMIDMap data structure
    load_way_id();
        
    // Load the intersectionStreetSegments and Load the intersectionStreetNames    
    load_street_segments_and_street_names_of_intersection(numIntersections);
    
    // Load the time spent travelling along every street segment
    load_time_taken(numStreetSegments);
    
    // Load the adjacentIntersections
    load_adjacent_intersections(numIntersections);
    
    // Populate the streetNamesDatabase
    std::string currentStreetName;
    for(StreetIndex streetID = 0; streetID < numStreets; streetID++){
        currentStreetName = getStreetName(streetID);
        removeSpacesAndCaptialsInString(currentStreetName);
        streetNamesDatabase.insert(std::make_pair(currentStreetName, streetID));
    }
    
    return true; // Loading was successful
}

void close_map() {
    // Clean up all the data structures that were loaded in load_map())
    streetSegmentsOfStreets.clear();
    intersectionsAlongStreets_v.clear();
    intersectionsAlongStreets_uos.clear();
    intersectionStreetSegments.clear();
    nodeIDToOSMIDMap.clear();
    total_travel_time.clear();
    wayIDToOSMIDMap.clear();
    intersectionStreetNames.clear();
    adjacentIntersections.clear();
    streetNamesDatabase.clear();
    
    // Close the StreetsDatabaseAPI and OSMDatabaseAPI
    closeStreetDatabase();
    closeOSMDatabase();
}

// Returns the rectangular projection of a Latitude/Longitude pair
std::pair<double, double> getEquirectangularProjection(LatLon point, double avgLatitude) {
    return std::make_pair(point.lon() * DEGREE_TO_RADIAN * std::cos(avgLatitude * DEGREE_TO_RADIAN), point.lat() * DEGREE_TO_RADIAN);
}

double find_distance_between_two_points(std::pair<LatLon, LatLon> points) {// uses distance formula to find distance between two points after converting to radians
    std::pair<double,double>x;
    std::pair<double,double>y;
    double distance;
    x.first=points.first.lon() * DEGREE_TO_RADIAN * std::cos((points.first.lat() * DEGREE_TO_RADIAN + points.second.lat() * DEGREE_TO_RADIAN)/2);
    y.first= points.first.lat() * DEGREE_TO_RADIAN;
    x.second=points.second.lon() * DEGREE_TO_RADIAN * std::cos((points.first.lat() * DEGREE_TO_RADIAN + points.second.lat() * DEGREE_TO_RADIAN)/2);    
    y.second= points.second.lat() * DEGREE_TO_RADIAN;
    
    distance = EARTH_RADIUS_METERS * std::sqrt(std::pow(y.second-y.first,2) + std::pow(x.second-x.first,2));
    return distance;
}
double find_street_segment_length(int street_segment_id) {

    double street_seg_length;
    InfoStreetSegment street_info = getInfoStreetSegment(street_segment_id);//get info of street whos legnth to find
    LatLon start_position = getIntersectionPosition(street_info.from);//find the from and to points to go to
    LatLon end_position = getIntersectionPosition(street_info.to);
    
    if(street_info.curvePointCount ==0) {// if there is no curve between find and to, find direct distance
        street_seg_length = (find_distance_between_two_points(std::make_pair(start_position,end_position)));
    }
    
    else {// if there are curves/turns to take, we find the sum of all straight segments and sum them all
        int j=0;
        double total = 0;
        LatLon previous_curve_position = start_position;
        LatLon current_curve_position; 
        
        for (j=0; j<=street_info.curvePointCount;j++){
            if(j== street_info.curvePointCount){
                current_curve_position=end_position;
            }
            else{
                current_curve_position =getStreetSegmentCurvePoint(j,street_segment_id);
            }
            
            total = find_distance_between_two_points(std::make_pair(previous_curve_position,current_curve_position))+ total;
            previous_curve_position=current_curve_position;
        }
        
        street_seg_length=total;       
    }
    
    return street_seg_length;
}

double find_street_segment_travel_time(int street_segment_id) {
    return total_travel_time[street_segment_id];// uses load map for the calculations and finds time using speed = dist/time;
}


int find_closest_intersection(LatLon my_position) {
    int intersection = 0;
    LatLon intersection2 = getIntersectionPosition(0);
    double d2 = pow((intersection2.lat() - my_position.lat()), 2) + pow((intersection2.lon() - my_position.lon()), 2) + 1;
    double cosine = cos(my_position.lat() * DEGREE_TO_RADIAN);
    double x2 = my_position.lon() * cosine;
    // Calculate the distance between every intersection point and the given point
    // Compare these distances and find the minimum one, then return the id of this intersection
    for (int n = 0; n < getNumIntersections(); n++) {
        LatLon intersection1 = getIntersectionPosition(n);
        double x1 = intersection1.lon() * cosine;
        double a2 = pow((intersection1.lat() - my_position.lat()), 2) + pow((x2 - x1), 2);
        if (a2 < d2) {// if we find another minimum one we replace it
            d2 = a2;
            intersection = n;
        }
    }
    return intersection;
}


//Returns the street segments for the given intersection
//Create an empty vector, loop through all streets in intersections
std::vector<int> find_street_segments_of_intersection(int intersection_id){
    return intersectionStreetSegments[intersection_id];
}

//Returns the street names at the given intersection (includes duplicate street 
//names in returned vector)
//This may change the global variable!
std::vector<std::string> find_street_names_of_intersection(int intersection_id) {
    return intersectionStreetNames[intersection_id];
}

//Compare each intersection list of roads with each other and check for one-way roads.
bool are_directly_connected(std::pair<int, int> intersectionIDs) {
    //Handle corner case of intersection connected to itself, same intersection passed in
    if (intersectionIDs.first == intersectionIDs.second)
        return true;
    
    std::vector<int> intersectionStreetSegments1 = find_street_segments_of_intersection(intersectionIDs.first);
    StreetSegmentIndex currentRoad;
    InfoStreetSegment currentRoadInfo;
    
    //Loop through all street segments of first street, compare to second id passed in
    for(int index1=0; index1 < intersectionStreetSegments1.size(); index1++){
        currentRoad = getIntersectionStreetSegment(intersectionIDs.first, index1); //segment number is index1 (from road1->road2))
        currentRoadInfo = getInfoStreetSegment(currentRoad);
        if (currentRoadInfo.oneWay == true){
            //Only check one case where points can be connected
            if(currentRoadInfo.from == intersectionIDs.first && currentRoadInfo.to == intersectionIDs.second){
                return true;
            }
        } else if (currentRoadInfo.oneWay == false){
            //Handle both cases where road is 2 way but still see if it connects either way
            if(currentRoadInfo.from == intersectionIDs.first && currentRoadInfo.to == intersectionIDs.second){
                return true;
            } else if (currentRoadInfo.to == intersectionIDs.first && currentRoadInfo.from == intersectionIDs.second){
                return true;
            }
        }
    }
    return false;
}

std::vector<int> find_adjacent_intersections(int intersection_id) {
    return adjacentIntersections[intersection_id];
}

std::vector<int> find_street_segments_of_street(int street_id) {  
    // This data is already organized into the streetSegmentsOfStreets vector
    // Simply fetch the data corresponding to the input street_id
    return streetSegmentsOfStreets[street_id];
}

std::vector<int> find_intersections_of_street(int street_id) {
    // This data is already organized into the intersectionsAlongStreets_v vector
    // Simply fetch the data corresponding to the input street_id
    return intersectionsAlongStreets_v[street_id];
}

std::vector<int> find_intersections_of_two_streets(std::pair<int, int> street_ids) {
    // Fetch an unordered_set consisting of all the intersections that lie on
    // each of the two passed-in streets
    //
    // Note: first one is fetched as a vector since it will be traversed sequentially
    std::vector<IntersectionIndex> intersectionsOfFirstStreet = intersectionsAlongStreets_v[street_ids.first];
    // Note: second one is fetched as an unordered_set since we will be searching
    // it by value
    std::unordered_set<IntersectionIndex> intersectionsOfSecondStreet = intersectionsAlongStreets_uos[street_ids.second];
    
    // List of common intersections that will be populated below
    std::vector<IntersectionIndex> commonIntersectionsOfStreets;
    
    // Iterate through all intersections in the first street, and check if they
    // also exist in the second street. 
    for (int n_intersection = 0; n_intersection < intersectionsOfFirstStreet.size(); n_intersection++) {
        // If the n'th intersection of the first street exists in the second, 
        // add it to the list of common intersections
        if (intersectionsOfSecondStreet.find(intersectionsOfFirstStreet[n_intersection]) != intersectionsOfSecondStreet.end()) {
            commonIntersectionsOfStreets.push_back(intersectionsOfFirstStreet[n_intersection]);
        }
    }

    // Return the (now populated) list of common intersections
    return commonIntersectionsOfStreets;
}

// Converts the passed-in string to all lowercase and removes any spaces
void removeSpacesAndCaptialsInString(std::string& str) {
    int spacesRemoved = 0;
    
    for (int i = 0; i < str.length(); i++) {
        if (str[i] == ' ') {
            spacesRemoved++;
        }
        else {
            str[i - spacesRemoved] = tolower(str[i]);
        }
    }
    
    str.resize(str.size() - spacesRemoved);
}

std::vector<int> find_street_ids_from_partial_street_name(std::string street_prefix){
    std::vector<int> possibleStreets;
      
    removeSpacesAndCaptialsInString(street_prefix);

    if(street_prefix.length()==0 || street_prefix == ""){// case if nothing is entered, return empty vector
        return possibleStreets;
    }
    
    removeSpacesAndCaptialsInString(street_prefix);// removes spaces and convert to lowercase

    std::string incremented = street_prefix;

    incremented[incremented.size() - 1]++;//use a map to arrange names alphabetically and search the name
  //use lower bound to find such iterators with match the first and last key of the map matching the street prefix argument
    auto matchLower = streetNamesDatabase.lower_bound(street_prefix);
    auto matchUpper = streetNamesDatabase.lower_bound(incremented);

    for(auto match = matchLower; match != matchUpper; match++){
        possibleStreets.push_back(match -> second);
    }
    
    return possibleStreets;
}

// Given a vector of LatLon points, returns the average latitude of all points
double findAverageLatitude(std::vector<LatLon> points) {
    double total = 0;
    
    for (int i = 0; i < points.size(); i++) {
        total += points[i].lat();
    }
    
    return total / points.size();
}

double find_feature_area(int feature_id) {
    int numPointsInFeature = getFeaturePointCount(feature_id);
    
    // Create a list of all points that define the feature
    std::vector<LatLon> pointsOfFeatureLatLon;
    pointsOfFeatureLatLon.reserve(numPointsInFeature);
    
    for (int i = 0; i < numPointsInFeature; i++) {
        LatLon featurePointLatLon = getFeaturePoint(i, feature_id);   
        pointsOfFeatureLatLon.push_back(featurePointLatLon);
    }
    
    // Check if the polygon is actually closed (start point = end point)
    if ((pointsOfFeatureLatLon[numPointsInFeature - 1].lat() != pointsOfFeatureLatLon[0].lat())
            || (pointsOfFeatureLatLon[numPointsInFeature - 1].lon() != pointsOfFeatureLatLon[0].lon())) {
        return 0;
    }
    
    // Convert the latitude/longitude points in the above vector into 
    // coordinates on an equirectangular projection and store in a new vector
    std::vector<std::pair<double, double>> pointsOfFeatureXY;
    pointsOfFeatureXY.reserve(numPointsInFeature);
    
    double averageLatitude = findAverageLatitude(pointsOfFeatureLatLon);
    
    for (int i = 0; i < numPointsInFeature; i++) {
        pointsOfFeatureXY.push_back(getEquirectangularProjection(pointsOfFeatureLatLon[i], averageLatitude));
    }
    
    // Calculate the feature area (in the equirectangular projection coordinate 
    // system) using the Shoelace Algorithm
    double featureArea = 0;
    
    for (int i = 0; i < numPointsInFeature - 1; i++) {
        featureArea += (pointsOfFeatureXY[i].first * pointsOfFeatureXY[i + 1].second) - (pointsOfFeatureXY[i + 1].first * pointsOfFeatureXY[i].second);
    }
    
    featureArea += (pointsOfFeatureXY[numPointsInFeature - 1].first * pointsOfFeatureXY[0].second) - (pointsOfFeatureXY[0].first * pointsOfFeatureXY[numPointsInFeature - 1].second);
    featureArea = 0.5 * std::abs(featureArea);
    
    // Unit conversion - lengths in the equirectangular projection scale by
    // EARTH_RADIUS, so areas scale by (EARTH_RADIUS ^ 2)
    featureArea = EARTH_RADIUS_METERS * EARTH_RADIUS_METERS * featureArea;
    
    return featureArea;
}

double find_way_length(OSMID way_OSMID) {
    // Retrieve a pointer to the OSMWay object specified by way_id,
    // using the global wayIDToOSMIDMap data structure
    int way_id = wayIDToOSMIDMap[way_OSMID];
    const OSMWay* specifiedWay = getWayByIndex(way_id);
    
    // Running total storing the length of the way
    double totalLengthOfWay = 0;
    
    // Iterate through all the nodes on the specified OSMWay, calculate the
    // distance between adjacent nodes, and add it to the running total above
    const std::vector<OSMID>& nodesOfWay = getWayMembers(specifiedWay);
    
    for (int n_node = 0; n_node < nodesOfWay.size() - 1; n_node++) {
        // Get the OSMNode object for the n'th node on the specifiedWay
        int idNode_n = nodeIDToOSMIDMap[nodesOfWay[n_node]];
        const OSMNode* node_n = getNodeByIndex(idNode_n);
        
        // Get the OSMNode object for the (n+1)'th node on the specifiedWay
        int idNode_n_next = nodeIDToOSMIDMap[nodesOfWay[n_node + 1]];
        const OSMNode* node_n_next = getNodeByIndex(idNode_n_next);
        
        // Create a pair storing the LatLon coordinates of the two adjacent 
        // nodes on the way
        std::pair<LatLon, LatLon> coordsOfAdjacentNodesOnWay;
        coordsOfAdjacentNodesOnWay.first = getNodeCoords(node_n);
        coordsOfAdjacentNodesOnWay.second = getNodeCoords(node_n_next);
        
        // Calculate the length of the segment of the way (represented by the
        // two adjacent nodes above) using the find_distance_between_two_points()
        // function from above, and add it to the running total
        totalLengthOfWay += find_distance_between_two_points(coordsOfAdjacentNodesOnWay);
    }
    
    // By the end of the above loop, the total length has been found. Return it
    return totalLengthOfWay;
}
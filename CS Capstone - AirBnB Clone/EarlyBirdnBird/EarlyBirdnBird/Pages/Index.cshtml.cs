using DataAccess;
using Microsoft.AspNetCore.Mvc;
using Microsoft.AspNetCore.Mvc.RazorPages;
using Infrastructure.Models;
using Microsoft.EntityFrameworkCore;
using System;
using System.Linq;
using System.Threading.Tasks;
using Microsoft.AspNetCore.Mvc.Rendering;

namespace EarlyBirdnBird.Pages
{
    public class IndexModel : PageModel
    {
        private readonly ILogger<IndexModel> _logger;
        private readonly ApplicationDbContext _dbContext;

        private int PageSize = 6;
        public int CurrentPage { get; set; } = 1;
        public int TotalPages { get; set; }

        public IndexModel(ApplicationDbContext db)
        {
            _dbContext = db;
        }
        public IList<Rental> Rentals { get; set; }


        //Search Variables
        [BindProperty(SupportsGet = true)]
        public string? SearchString { get; set; }

        public SelectList? Bathrooms { get; set; }

        [BindProperty(SupportsGet = true)]
        public string? BathroomSearch { get; set; }

        [BindProperty(SupportsGet = true)]
        public string? BedSearch { get; set; }

        public SelectList? Occupants { get; set; }

        [BindProperty(SupportsGet = true)]
        public string? OccupantSearch { get; set; }

        [BindProperty(SupportsGet = true)]
        public double? PriceSearch { get; set; }

        [BindProperty(SupportsGet = true)]
        public DateTime? StartDate { get; set; }

        [BindProperty(SupportsGet = true)]
        public DateTime? EndDate { get; set; }
        
        public SelectList? Amenities { get; set; }

        [BindProperty(SupportsGet = true)]
        public string AmenitySearch {  get; set; }

        public SelectList? Locations { get; set; }

        [BindProperty(SupportsGet = true)]
        public string LocationSearch { get; set; }

        public SelectList? State { get; set; }

        [BindProperty(SupportsGet = true)]
        public string StateSearch { get; set; }

        public SelectList? Beds { get; set; }


        [BindProperty(SupportsGet = true)]
        public string? bedSearch { get; set; }

        public async Task OnGetAsync(int currentPage = 1)
        {
            CurrentPage = currentPage;

            //initial rental query, which will be slowly widdled down based on filters
            var rentals = from r in _dbContext.Rentals where r.Published == true select r;

            //Bathroom Query, gets the Bathrooms and orders them
            var bathroomQuery = from r in _dbContext.Rentals
                                orderby r.Bathrooms
                                select r.Bathrooms;

            Bathrooms = new SelectList(await bathroomQuery.Distinct().ToListAsync());

            //Occupant Query, gets the Occupant Max and orders them
            var OccuQuery = from r in _dbContext.Rentals
                            orderby r.MaxOcupants
                            select r.MaxOcupants;

            Occupants = new SelectList(await OccuQuery.Distinct().ToListAsync());

            //Occupant Query, gets the Occupant Max and orders them
            var bedQuery = from r in _dbContext.Rentals
                           orderby r.Beds
                           select r.Beds;

            Beds = new SelectList(await bedQuery.Distinct().ToListAsync());

            //Amenities Query, gets the Amenities from the RentalAmenitites Table, connects it all, and orders by name
            IQueryable<string> amenitiyQuery = from a in _dbContext.RentalAmenities
                                               join r in _dbContext.Rentals on a.RentalId equals r.Id
                                               join i in _dbContext.Amenities on a.AmenityId equals i.Id
                                               orderby i.Name
                                               select i.Name;
            Amenities = new SelectList(await amenitiyQuery.Distinct().ToListAsync());

            //Location Query, gets the Cities and orders them
            IQueryable<string> locationQuery = from r in _dbContext.Rentals
                                               orderby r.City
                                               select r.City;
            Locations = new SelectList(await locationQuery.Distinct().ToListAsync());

            //State Query, gets the States and orders them
            IQueryable<string> stateQuery = from r in _dbContext.Rentals
                                            orderby r.State
                                            select r.State;

            State = new SelectList(await stateQuery.Distinct().ToListAsync());

            //Actual filtering beyond here

            if (!string.IsNullOrEmpty(SearchString))
            {
                rentals = rentals.Where(s => s.Name.Contains(SearchString));
            }

            if (!string.IsNullOrEmpty(BedSearch))
            {
                int bedNum = Convert.ToInt32(BedSearch);
                rentals = rentals.Where(s => s.Beds == bedNum);
            }

            if (!string.IsNullOrEmpty(BathroomSearch))
            {
                int bathroomNum = Convert.ToInt32(BathroomSearch);
                rentals = rentals.Where(s => s.Bathrooms == bathroomNum);
            }

            if (!string.IsNullOrEmpty(OccupantSearch))
            {
                int peopleNum = Convert.ToInt32(OccupantSearch);
                rentals = rentals.Where(s => s.MaxOcupants == peopleNum);
            }

            if (!string.IsNullOrEmpty(AmenitySearch))
            {
                var amenQuery = from a in _dbContext.RentalAmenities
                                join r in _dbContext.Rentals on a.RentalId equals r.Id
                                join i in _dbContext.Amenities on a.AmenityId equals i.Id
                                where i.Name == AmenitySearch.ToString()
                                select r.Id;

                rentals = rentals.Where(s => amenQuery.Contains(s.Id));
            }

            if (!string.IsNullOrEmpty(LocationSearch))
            {
                rentals = rentals.Where(s => s.City.Contains(LocationSearch));
            }

            if (!string.IsNullOrEmpty(StateSearch))
            {
                rentals = rentals.Where(s => s.State.Contains(StateSearch));
            }

            //if (PriceSearch.HasValue) 
            //{
            //    int price = Convert.ToInt32(BathroomSearch);
            //    rentals = rentals.Where(s => s.Prices.Any(p => p.PricePerNight == price));
            //}

            //checks if there are any reservations for the selected dates. 
            //if there are no reservations for the selected dates, then it will return a rental.
            if (StartDate.HasValue && EndDate.HasValue)
            {
                rentals = rentals.Where(s => !_dbContext.Reservation.Any(r =>
                    r.Rental == s.Id &&
                    ((r.BookedStart <= StartDate && r.BookedEnd >= StartDate) ||
                    (r.BookedStart <= EndDate && r.BookedEnd >= EndDate)) &&
                    (r.Status == 1 || r.Status == 2)));
            }

            int totalRentals = await rentals.CountAsync();
            TotalPages = (int)Math.Ceiling(totalRentals / (double)PageSize);

            Rentals = await rentals
                .Skip((CurrentPage - 1) * PageSize)
                .Take(PageSize)
                .Include(r => r.Images)
                .Include(r => r.Prices)
                .Include(r => r.OwnerName)
            .ToListAsync();
            //Bathrooms = new SelectList(await bathroomQuery.Distinct().ToListAsync());
            //Occupants = new SelectList(await OccupantQuery.Distinct().ToListAsync());
        }
    }
}
using DataAccess;
using Infrastructure.Models;
using Microsoft.AspNetCore.Mvc.RazorPages;
using Microsoft.EntityFrameworkCore;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace EarlyBirdnBird.Pages
{
    public class HostSummaryModel : PageModel
    {
        private readonly ApplicationDbContext _dbContext;

        public HostSummaryModel(ApplicationDbContext dbContext)
        {
            _dbContext = dbContext;
        }

        public IList<Reservation> Reservations { get; set; }
        public double TotalEarnings { get; set; }
        public int TotalBookings { get; set; }
        public double AverageEarningsPerBooking { get; set; }

        public async Task OnGetAsync()
        {
            DateTime today = DateTime.Now;
//ADJUSTED FOR DB CHANGES
            Reservations = await _dbContext.Reservation
                .Include(r => r.AppUser)
                .Where(r => r.BookedStart <= today && today <= r.BookedEnd)
                .ToListAsync();

            TotalEarnings = Reservations.Sum(r => r.GrandTotal);
            TotalBookings = Reservations.Count;
            AverageEarningsPerBooking = TotalBookings > 0 ? TotalEarnings / TotalBookings : 0;
        }
    }
}